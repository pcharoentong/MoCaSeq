#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
	extract_data;
	file_has_extension
} from "./modules/input"

include {
	parse_stub_json
} from "./modules/stub"

include {
	PREPARE_GENOME;
	GENOME_ANNOTATION
} from "./modules/local/subworkflow/genome"

include {
	MAP as MAPPER
	REMAP as REMAPPER
} from "./modules/local/subworkflow/remap"

include {
	MANTA
} from "./modules/local/subworkflow/manta"

include {
	STRELKA
} from "./modules/local/subworkflow/strelka"

include {
	MUTECT
} from "./modules/local/subworkflow/mutect"

include {
	DELLY
} from "./modules/local/subworkflow/delly"

include {
	CNV_KIT;
	CNV_KIT_SELECT_SAMPLE;
	CNV_KIT_TARGET_BED;
	CNV_KIT_COVERAGE;
	CNV_KIT_FIX;
	CNV_KIT_SEGMENT;
	CNV_KIT_PON
} from "./modules/local/subworkflow/cnv-kit"

include {
	HMM_COPY
} from "./modules/local/subworkflow/hmm-copy"

include {
	LOH
} from "./modules/local/subworkflow/loh"

include {
	MSI_SENSOR
} from "./modules/local/subworkflow/msi-sensor"

include {
	BUBBLE_TREE
} from "./modules/local/subworkflow/bubble-tree"

include {
	JABBA
} from "./modules/local/subworkflow/jabba"

include {
	IGV_TRACK_READ;
	IGV_TRACK_CNR as IGV_TRACK_CNR_cnv_kit;
	IGV_TRACK_CNS;
	IGV_TRACK_CNS as IGV_TRACK_CNS_cnv_kit;
	IGV_TRACK_VCF_SV as IGV_TRACK_VCF_SV_jabba;
	IGV_TRACK_VCF_SV as IGV_TRACK_VCF_SV_manta
} from "./modules/local/subworkflow/igv-track"

include {
	COHORT_QC as COHORT_QC_human;
	COHORT_QC as COHORT_QC_mouse
} from "./modules/local/subworkflow/qc"

tsv_path = null
pon_bed_path = null
pon_tsv_path = null

ch_input_sample = Channel.empty ()


// check if we have valid --input
if (params.input == null && params.pon_tsv == null && params.qc_dir == null) {
	  exit 1, "[MoCaSeq] error: --input or --pon_tsv or --qc_dir were not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
if (params.input && (file_has_extension (params.input, "tsv"))) tsv_path = params.input
if (params.pon_bed && (file_has_extension (params.pon_bed, "bed"))) pon_bed_path = params.pon_bed
if (params.pon_tsv && (file_has_extension (params.pon_tsv, "tsv"))) pon_tsv_path = params.pon_tsv

if (tsv_path) {

	tsv_file = file (tsv_path)
	if (tsv_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[MoCaSeq] error: input TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.input}"
	ch_input_sample = extract_data (tsv_path)

}
else if (pon_tsv_path)
{
	tsv_file = file (pon_tsv_path)
	if (tsv_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one TSV file per run."
	if (!tsv_file.exists ()) exit 1, "[MoCaSeq] error: pon_tsv TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.pon_tsv}"
}
else if (params.qc_dir)
{
	qc_dir_file = file (params.qc_dir)
	if (qc_dir_file instanceof List) exit 1, "[MoCaSeq] error: can only accept one QC dir per run."
	if ( ! ( qc_dir_file.exists () && qc_dir_file.isDirectory () ) ) exit 1, "[MoCaSeq] error: qc_dir directory could not be found or was not a directory. You gave the path: ${params.qc_dir}"
}
else exit 1, "[MoCaSeq] error: --input or --pon_tsv file(s) or --qc_dir not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."

// Optionally load json map to control the behaviour of stubs (cp vs touch)
if (params.stub_json && ( file_has_extension (params.stub_json, "js") || file_has_extension (params.stub_json, "json") ) ) {
	params.stub_json_map = parse_stub_json (params.stub_json)
}

ch_input_branched = ch_input_sample.branch {
	bam: it["normalBAM"] != null || it["tumorBAM"] != null //These are all BAMs
	map: ( it["normalR1"] != null && it["normalR1"].toString().endsWith (".fastq.gz") ) || ( it["tumorR1"] != null && it["tumorR1"].toString().endsWith (".fastq.gz") ) //Path.endsWith tries to match entire final segment
	remap: ( it["normalR1"] != null && it["normalR1"].toString().endsWith (".bam") ) || ( it["tumorR1"] != null && it["tumorR1"].toString().endsWith (".bam") ) //Path.endsWith tries to match entire final segment
}

ch_input_branched_bam_branched = ch_input_branched.bam.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_map_branched = ch_input_branched.map.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wgs: it["organism"] == "mouse" && it["seqType"] == "wgs"
	other: true
}

ch_input_branched_remap_branched = ch_input_branched.remap.branch {
	human_wgs: it["organism"] == "human" && it["seqType"] == "wgs"
	mouse_wgs: it["organism"] == "mouse" && it["seqType"] == "wgs"
	mouse_wex: it["organism"] == "mouse" && it["seqType"] == "wex"
	other: true
}

ch_input_branched_bam_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input bam:\n${it}" }
ch_input_branched_map_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input map:\n${it}" }
ch_input_branched_remap_branched.other.view { "[MoCaSeq] error: Failed to find matching workflow (organism & seqType) for input remap:\n${it}" }

workflow HUMAN_WGS
{
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	ch_bam = ch_input_branched_bam_branched.human_wgs

	MANTA (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam)
	STRELKA (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
	MUTECT (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	DELLY (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_bam)
	HMM_COPY (params.genome_build.human, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, MUTECT.out.result)
	MSI_SENSOR (params.genome_build.human, GENOME_ANNOTATION.out.micro_satellite, ch_bam)

	if ( params.pon_dir == null )
	{
		CNV_KIT (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index_flat, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gencode_genes_bed, ch_bam)
		BUBBLE_TREE (params.genome_build.human, PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.call, LOH.out.result)
		JABBA (params.genome_build.human, PREPARE_GENOME.out.chrom_names, MANTA.out.vcf, HMM_COPY.out.cnr, HMM_COPY.out.call, BUBBLE_TREE.out.result)
	}
	else
	{
		ch_bam_tumor = ch_bam.filter { it["type"] == "Tumor" }.map { tuple (it, "Tumor", it["tumorBAM"], it["tumorBAI"] ) }
		ch_target_bed = Channel.of ([
			file ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.target.bed", glob: false),
			file ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.resolution.json", glob: false)
		]).map {
			def jsonSlurper = new groovy.json.JsonSlurper ()
			def data = jsonSlurper.parse (it[1])
			tuple ( it[0], data["target_avg_size"], data["wgs_depth"] )
		}.first ()

		CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_target_bed, ch_bam_tumor)
		CNV_KIT_FIX (params.genome_build.human, Channel.fromPath ("${params.pon_dir}/${params.genome_build.human}_PON/${params.genome_build.human}.reference.cnn").first (), CNV_KIT_COVERAGE.out.result)
		CNV_KIT_SEGMENT (params.genome_build.human, CNV_KIT_FIX.out.cnr)
		BUBBLE_TREE (params.genome_build.human, PREPARE_GENOME.out.chrom_names_auto, CNV_KIT_SEGMENT.out.call, LOH.out.result)
		JABBA (params.genome_build.human, PREPARE_GENOME.out.chrom_names, MANTA.out.vcf, CNV_KIT_FIX.out.cnr, CNV_KIT_SEGMENT.out.call, BUBBLE_TREE.out.result)

		if ( params.track_cn )
		{
			IGV_TRACK_CNR_cnv_kit (params.genome_build.human,  PREPARE_GENOME.out.interval_bed, CNV_KIT_FIX.out.cnr)
			IGV_TRACK_CNS_cnv_kit (params.genome_build.human, CNV_KIT_SEGMENT.out.cns)
			IGV_TRACK_CNS_hmm_copy (params.genome_build.human, HMM_COPY.out.cnr)
		}
		if ( params.track_sv )
		{
			IGV_TRACK_VCF_SV_jabba (params.genome_build.human, JABBA.out.vcf.map { [ it[0], "JaBbA", it[1] ] })
			IGV_TRACK_VCF_SV_manta (params.genome_build.human, MANTA.out.vcf.map { [ it[0], "Manta", it[1] ] })
		}
	}

	if ( params.track_read )
	{
		IGV_TRACK_READ (params.genome_build.human, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, ch_bam)
	}
	if ( params.track_cn )
	{
		IGV_TRACK_CNS (params.genome_build.human, CNV_KIT.out.cns)
	}
}

workflow MOUSE_WEX
{
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	ch_bam = ch_input_branched_bam_branched.mouse_wex

	MANTA (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam)
	STRELKA (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam, MANTA.out.indel)
	MUTECT (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out._chrom_n, ch_bam)
	HMM_COPY (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, GENOME_ANNOTATION.out.gc_wig, GENOME_ANNOTATION.out.map_wig, ch_bam)
	LOH (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fasta_index, PREPARE_GENOME.out.chrom_names, PREPARE_GENOME.out.interval_bed, MUTECT.out.result)
	MSI_SENSOR (params.genome_build.mouse, GENOME_ANNOTATION.out.micro_satellite, ch_bam)
	BUBBLE_TREE (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names_auto, HMM_COPY.out.call, LOH.out.result)
}

workflow HUMAN_MAP {
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	MAPPER (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_map_branched.human_wgs)
	REMAPPER (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.human_wgs)
}

workflow MOUSE_MAP {
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	MAPPER (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_map_branched.mouse_wgs)
	REMAPPER (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.dir, GENOME_ANNOTATION.out.common_vcf, ch_input_branched_remap_branched.mouse_wgs)
}

workflow HUMAN_PON {
	main:
	PREPARE_GENOME (params.genome_build.human)
	GENOME_ANNOTATION (params.genome_build.human)

	if ( pon_bed_path == null )
	{
		// Generate target regions for CNVKit
		ch_bam_normal = ch_input_branched_bam_branched.human_wgs.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		if ( params.pon_sample == null )
		{
			// First we need to pick a median coverage bam
			ch_normal_size_lines = ch_bam_normal.map { [ it[0]["sampleName"], params.genome_build.human, it[0]["normalBAM"], it[0]["normalBAI"], it[0]["normalBAM"].size () ].join ("\t") }
			ch_normal_size_tsv = Channel.of ( ["sample", "genome_build", "bam_path", "bai_path", "bam_size"].join ("\t") )
				.concat (ch_normal_size_lines)
				.collectFile (name: "${params.genome_build.human}.normal_sizes.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

			if ( params.pon_dry )
			{
				CNV_KIT_SELECT_SAMPLE (params.genome_build.human, ch_normal_size_tsv)
				ch_bam_normal_sample = ch_bam_normal.map { [ it[0]["sampleName"], it] }.join (CNV_KIT_SELECT_SAMPLE.out.result).map { it[1] }
				CNV_KIT_TARGET_BED (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, ch_bam_normal_sample)
				CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, CNV_KIT_TARGET_BED.out.result, ch_bam_normal)

				ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.human, it[2], it[3]].join ("\t") }
				ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
					.concat (ch_normal_coverage_lines)
					.collectFile (name: "${params.genome_build.human}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

				CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_normal_coverage_tsv)
			}
		}
		else
		{
			ch_bam_normal_sample = ch_bam_normal.first { it[0]["sampleName"] == params.pon_sample }
			CNV_KIT_TARGET_BED (params.genome_build.human, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.chrom_names, ch_bam_normal_sample)
		}
	}
	else
	{
		if ( pon_tsv_path == null )
		{
			if ( params.pon_resolution == null ) exit 1, "[MoCaSeq] error: You must also supply --pon_resolution when you give --pon_tsv (to ensure consistent names for the coverage .cnn files)"
			ch_bam_normal = ch_input_branched_bam_branched.human_wgs.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }
			ch_target_bed = Channel.value ( [ file (pon_bed_path, glob: false), params.pon_resolution as long, 0 ] )

			CNV_KIT_COVERAGE (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_target_bed, ch_bam_normal)

			ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.human, it[2], it[3]].join ("\t") }
			ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
				.concat (ch_normal_coverage_lines)
				.collectFile (name: "${params.genome_build.human}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.human}_PON")

			if ( params.pon_dry )
			{
				CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, ch_normal_coverage_tsv)
			}
		}
		else
		{
			CNV_KIT_PON (params.genome_build.human, PREPARE_GENOME.out.fasta, Channel.fromPath (pon_tsv_path))
		}
	}
}

workflow MOUSE_PON {
	main:
	PREPARE_GENOME (params.genome_build.mouse)
	GENOME_ANNOTATION (params.genome_build.mouse)

	if ( pon_tsv_path == null )
	{
		ch_bam_normal = ch_input_branched_bam_branched.mouse_wex.filter { it["type"] == "Normal" }.map { tuple (it, "Normal", it["normalBAM"], it["normalBAI"] ) }

		CNV_KIT_COVERAGE (params.genome_build.mouse, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.interval_bed, ch_bam_normal)

		ch_normal_coverage_lines = CNV_KIT_COVERAGE.out.result.map { [it[0]["sampleName"], params.genome_build.mouse, it[2], it[3]].join ("\t") }
		ch_normal_coverage_tsv = Channel.of ( ["sample", "genome_build", "resolution", "normal_cov"].join ("\t")  )
			.concat (ch_normal_coverage_lines)
			.collectFile (name: "${params.genome_build.mouse}.normal_coverage_file_paths.tsv", newLine: true, sort: false, storeDir: "${params.output_base}/${params.genome_build.mouse}_PON")

		if ( params.pon_dry )
		{
			//CNV_KIT_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, ch_normal_coverage_tsv)
		}
	}
	else
	{
		//CNV_KIT_PON (params.genome_build.mouse, PREPARE_GENOME.out.chrom_names, GENOME_ANNOTATION.out.par_interval_bed, Channel.fromPath (pon_tsv_path))
	}
}

workflow {
	HUMAN_WGS ()
	MOUSE_WEX ()
}

// Run using -entry MAP
workflow MAP {
	HUMAN_MAP ()
	MOUSE_MAP ()
}

// Run using -entry PON
workflow PON {
	HUMAN_PON ()
	MOUSE_PON ()
}


// Run using -entry QC
workflow QC {
	COHORT_QC_human (params.genome_build.human, Channel.fromPath (params.qc_dir))
	COHORT_QC_mouse (params.genome_build.mouse, Channel.fromPath (params.qc_dir))
}

