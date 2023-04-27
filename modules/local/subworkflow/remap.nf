

include {
	sam_to_fastq_paired;
	fastqc_paired as fastqc_paired_extracted;
	fastqc_paired as fastqc_paired_trimmed;
	trim_paired;
	bwa_mem_paired;
	mark_duplicates;
	recalibrate
} from "../software/remap/main"

workflow MAP
{
	take:
		genome_build
		ch_bwa_index
		ch_fasta
		ch_common_vcf
		ch_data

	main:
	if (params.debug) { println "entered MAP subworfklow" }
		ch_data_branched = ch_data.map { it ->
			if ( it["type"] == "Normal" )
			{
				return (tuple (it, "Normal", it["normalR1"], it["normalR2"]))
			}
			if ( it["type"] == "Tumor" )
			{
				return (tuple (it, "Tumor", it["tumorR1"], it["tumorR2"]) )
			}
		}
		.dump (tag: 'remap input')
		.branch {
			paired: it[2].toString ().endsWith (".fastq.gz") && it[3].toString ().endsWith (".fastq.gz")
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching MAP workflow path for input:\n${it}" }

		fastqc_paired_extracted (genome_build, Channel.value ("MAP_input"), ch_data_branched.paired)
		trim_paired (fastqc_paired_extracted.out.result)
		fastqc_paired_trimmed (genome_build, Channel.value ("MAP_trimmed"), trim_paired.out.result)
		bwa_mem_paired (ch_fasta, fastqc_paired_trimmed.out.result)
		mark_duplicates (genome_build, bwa_mem_paired.out.result)
		recalibrate (genome_build, ch_fasta, ch_common_vcf, mark_duplicates.out.result)
		sample = recalibrate.out.result.map {
			it[0][it[0][[it[1].toLowerCase (),"BAM"].join ("")]] = it[2]
			it[0][it[0][[it[1].toLowerCase (),"BAI"].join ("")]] = it[3]
			it[0]
		}.reduce ( [:] ) { accumulator, item ->
			// Group by sample group
			if ( accumulator.containsKey (item["sampleGroup"]) )
			{
				accumulator[item["sampleGroup"]].add (item)
			}
			else
			{
				accumulator[item["sampleGroup"]] = [item]
			}
			accumulator
		}
		.flatMap ().flatMap { it ->
			// Within the sample group annotate Tumor samples with the new Normal bam
					def samples_by_type = it.value.inject ([:]) { accumulator, item ->
				if ( accumulator.containsKey (item["type"]) )
				{
					accumulator[item["type"]].add (item)
				}
				else
				{
					accumulator[item["type"]] = [item]
				}
				accumulator
			}
			println ("samples_by_type from map")
			println (samples_by_type)

			def result = []
			if ( samples_by_type.containsKey ("Normal") && samples_by_type.containsKey ("Tumor") )
			{
				result.addAll (samples_by_type["Normal"].collect { jt -> jt.putAll (["sampleGroupSize":it.value.size ()]);jt})
				result.addAll (samples_by_type["Tumor"].collect { jt -> jt.putAll (samples_by_type["Normal"][0].subMap (["normalBAM", "normalBAI"]));jt.put ("sampleGroupSize",it.value.size ());jt } )
			}
			else
			{
				samples_by_type.values ().each { jt -> result.addAll ( jt.collect { kt -> kt.put ("sampleGroupSize",it.value.size ());kt } ) }
			}
			result

		}.dump (tag: 'MAP output')

	emit:
		result = sample
}

workflow REMAP
{
	take:
		genome_build
		ch_bwa_index
		ch_fasta
		ch_common_vcf
		ch_data

	main:
	if (params.debug) { println "entered REMAP subworfklow" }
		ch_data_branched = ch_data.map { it ->
			if ( it["type"] == "Normal" )
			{
				return (tuple (it, "Normal", it["normalR1"], it["normalR2"]))
			}
			if ( it["type"] == "Tumor" )
			{
				return (tuple (it, "Tumor", it["tumorR1"], it["tumorR2"]) )
			}
		}
		.dump (tag: 'remap input')
		.branch {
			paired: it[2].toString ().endsWith (".bam") && it[2] == it[3]
			other: true
		}

		ch_data_branched.other.view { "[MoCaSeq] error: Failed to find matching REMAP workflow path for input:\n${it}" }
		
		ch_data_branched.paired.view { "args for sam_to_fastq_paired\n ${it[0]}, ${it[1]}, ${it[2]}" }
		// println "args for sam_to_fastq_paired" // ${it[0]}, ${it[1]}, ${it[2]}"
		sam_to_fastq_paired (ch_data_branched.paired.map { tuple (it[0], it[1], it[2] ) })
		sam_to_fastq_paired.out.result.view { "args for fastqc_paired_extracted\n ${it}" }
		println genome_build
		fastqc_paired_extracted (genome_build, Channel.value ("REMAP_extracted") , sam_to_fastq_paired.out.result)
		// println "args for bwa_mem_paired: "
		fastqc_paired_extracted.out.result.view {"args for bwa_mem_paired:\n ${it}"}
		ch_bwa_index.view ()
		bwa_mem_paired (ch_bwa_index, fastqc_paired_extracted.out.result)
		// println "args for mark_duplicates: "
		bwa_mem_paired.out.result.view {"args for mark_duplicates:\n ${it}"}
		mark_duplicates (genome_build, bwa_mem_paired.out.result)
		// println "args for recalibrate: "
		mark_duplicates.out.result.view {"args for recalibrate:\n ${it}"}
		recalibrate (genome_build, ch_fasta, ch_common_vcf, mark_duplicates.out.result)
		sample = recalibrate.out.result.map { [it[0]["sampleName"], it] }
			.groupTuple (size: 2)
			.map { it[1] }
			.map { it ->
				def bam_info = it.inject ([:]) { accumulator, item ->
					accumulator[[item[1].toLowerCase (),"BAM"].join ("")] = item[2]
					accumulator[[item[1].toLowerCase (),"BAI"].join ("")] = item[3]
					accumulator
				}
				def result = it[0][0].clone ()
				result.putAll (bam_info)
				result
			}

	emit:
		result = sample
}

