$HOSTNAME = ""
params.outdir = 'results'  

// global param
params.mate = "pair"
params.mate2 = "single"

// Process Parameters:
// Process Parameters for Assemble_pairs_assemble_pairs:
params.Assemble_pairs_assemble_pairs.method = "align" 
params.Assemble_pairs_assemble_pairs.coord = "illumina" 
params.Assemble_pairs_assemble_pairs.rc = "tail"
params.Assemble_pairs_assemble_pairs.head_fields_R1 = "UMI TRIM"
params.Assemble_pairs_assemble_pairs.head_fields_R2 = ""
params.Assemble_pairs_assemble_pairs.failed = "true"
params.Assemble_pairs_assemble_pairs.fasta = "false"
params.Assemble_pairs_assemble_pairs.nproc = "${params.nproc}"
params.Assemble_pairs_assemble_pairs.alpha = 0.00001
params.Assemble_pairs_assemble_pairs.maxerror = 0.3
params.Assemble_pairs_assemble_pairs.minlen = 8
params.Assemble_pairs_assemble_pairs.maxlen = 1000
params.Assemble_pairs_assemble_pairs.scanrev = "false"
params.Assemble_pairs_assemble_pairs.minident = 0.5
params.Assemble_pairs_assemble_pairs.evalue = 0.00001
params.Assemble_pairs_assemble_pairs.maxhits = 100
params.Assemble_pairs_assemble_pairs.fill = "false"
params.Assemble_pairs_assemble_pairs.aligner = "blastn"
params.Assemble_pairs_assemble_pairs.// align_exec =  ""
params.Assemble_pairs_assemble_pairs.// dbexec =  ""   
params.Assemble_pairs_assemble_pairs.gap = 0
params.Assemble_pairs_assemble_pairs.usearch_version = "11.0.667"
params.Assemble_pairs_assemble_pairs.assemble_reference =  ''

// Process Parameters for Assemble_pairs_parse_log_AP:
params.Assemble_pairs_parse_log_AP.field_to_parse = "ID REFID LENGTH OVERLAP GAP ERROR IDENTITY PVALUE EVALUE1 EVALUE2"

// Process Parameters for Assemble_pairs_ref_assemble_pairs:
params.Assemble_pairs_ref_assemble_pairs.method = "reference" 
params.Assemble_pairs_ref_assemble_pairs.coord = "illumina" 
params.Assemble_pairs_ref_assemble_pairs.rc = "tail"
params.Assemble_pairs_ref_assemble_pairs.head_fields_R1 = "UMI TRIM"
params.Assemble_pairs_ref_assemble_pairs.head_fields_R2 = ""
params.Assemble_pairs_ref_assemble_pairs.failed = "false"
params.Assemble_pairs_ref_assemble_pairs.fasta = "false"
params.Assemble_pairs_ref_assemble_pairs.nproc = "${params.nproc}"
params.Assemble_pairs_ref_assemble_pairs.alpha = 0.00001
params.Assemble_pairs_ref_assemble_pairs.maxerror = 0.3
params.Assemble_pairs_ref_assemble_pairs.minlen = 8
params.Assemble_pairs_ref_assemble_pairs.maxlen = 1000
params.Assemble_pairs_ref_assemble_pairs.scanrev = "false"
params.Assemble_pairs_ref_assemble_pairs.minident = 0.5
params.Assemble_pairs_ref_assemble_pairs.evalue = 0.00001
params.Assemble_pairs_ref_assemble_pairs.maxhits = 100
params.Assemble_pairs_ref_assemble_pairs.fill = "false"
params.Assemble_pairs_ref_assemble_pairs.aligner = "blastn"
params.Assemble_pairs_ref_assemble_pairs.// align_exec =  ""
params.Assemble_pairs_ref_assemble_pairs.// dbexec =  ""   
params.Assemble_pairs_ref_assemble_pairs.gap = 0
params.Assemble_pairs_ref_assemble_pairs.usearch_version = "11.0.667"
params.Assemble_pairs_ref_assemble_pairs.assemble_reference =  '/home/bcrlab/veredk/M24/IGH/ref.fasta'

// Process Parameters for Assemble_pairs_ref_parse_log_AP:
params.Assemble_pairs_ref_parse_log_AP.field_to_parse = "ID REFID LENGTH OVERLAP GAP ERROR IDENTITY PVALUE EVALUE1 EVALUE2"

// Process Parameters for Filter_Sequence_Quality_filter_seq_quality:
params.Filter_Sequence_Quality_filter_seq_quality.method = "quality"
params.Filter_Sequence_Quality_filter_seq_quality.nproc = "${params.nproc}"
params.Filter_Sequence_Quality_filter_seq_quality.q = "20"
params.Filter_Sequence_Quality_filter_seq_quality.n_length = "35"
params.Filter_Sequence_Quality_filter_seq_quality.n_missing = "10"


// Process Parameters for Filter_Sequence_Quality1_filter_seq_quality:
params.Filter_Sequence_Quality1_filter_seq_quality.method = "quality"
params.Filter_Sequence_Quality1_filter_seq_quality.nproc = "${params.nproc}"
params.Filter_Sequence_Quality1_filter_seq_quality.q = "20"
params.Filter_Sequence_Quality1_filter_seq_quality.n_length = "35"
params.Filter_Sequence_Quality1_filter_seq_quality.n_missing = "10"

if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mate2){params.mate2 = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_6_reads_g_11}
 } else {  
	g_6_reads_g_11 = Channel.empty()
 }

Channel.value(params.mate).into{g_7_mate_g_11;g_7_mate_g1_15;g_7_mate_g1_19;g_7_mate_g1_12;g_7_mate_g14_15;g_7_mate_g14_19;g_7_mate_g14_12}
Channel.value(params.mate2).into{g_10_mate_g2_7;g_10_mate_g2_5;g_10_mate_g2_0;g_10_mate_g19_7;g_10_mate_g19_5;g_10_mate_g19_0}


process unizp {

input:
 set val(name),file(reads) from g_6_reads_g_11
 val mate from g_7_mate_g_11

output:
 set val(name),file("*.fastq")  into g_11_reads0_g1_12

script:

readArray = reads.toString().split(' ')	
R1 = readArray.grep(~/.*R1.*/)[0]
R2 = readArray.grep(~/.*R2.*/)[0]

"""
case "$R1" in
*.gz | *.tgz ) 
        gunzip -c $R1 > R1.fastq
        ;;
*)
        cp $R1 ./R1.fastq
        echo "$R1 not gzipped"
        ;;
esac

case "$R2" in
*.gz | *.tgz ) 
        gunzip -c $R2 > R2.fastq
        ;;
*)
        cp $R2 ./R2.fastq
        echo "$R2 not gzipped"
        ;;
esac
"""
}


process Assemble_pairs_assemble_pairs {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_assemble-pass.f.*$/) "AP_align_read_pass/$filename"}
input:
 set val(name),file(reads) from g_11_reads0_g1_12
 val mate from g_7_mate_g1_12

output:
 set val(name),file("*_assemble-pass.f*")  into g1_12_reads0_g19_0
 set val(name),file("AP_*")  into g1_12_logFile1_g1_15
 set val(name),file("*_assemble-fail.f*") optional true  into g1_12_reads_failed2_g14_12
 set val(name),file("out*")  into g1_12_logFile3_g12_0

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi
	
	echo \$align_exec
	echo \$dbexec
	
	
	
	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc} >> out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

input:
 set val(name),file(log_file) from g1_12_logFile1_g1_15
 val mate from g_7_mate_g1_15

output:
 set val(name),file("*.tab")  into g1_15_logFile0_g1_19

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 set val(name),file(log_files) from g1_15_logFile0_g1_19
 val matee from g_7_mate_g1_19

output:
 file "*.rmd"  into g1_19_rMarkdown0_g1_22



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "assemble_pairs_report/$filename"}
input:
 file rmk from g1_19_rMarkdown0_g1_22

output:
 file "*.html"  into g1_22_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_ref_assemble_pairs {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_assemble-pass.f.*$/) "AP_ref_read_pass/$filename"}
input:
 set val(name),file(reads) from g1_12_reads_failed2_g14_12
 val mate from g_7_mate_g14_12

output:
 set val(name),file("*_assemble-pass.f*")  into g14_12_reads0_g2_0
 set val(name),file("AP_*")  into g14_12_logFile1_g14_15
 set val(name),file("*_assemble-fail.f*") optional true  into g14_12_reads_failed22
 set val(name),file("out*")  into g14_12_logFile33

script:
method = params.Assemble_pairs_ref_assemble_pairs.method
coord = params.Assemble_pairs_ref_assemble_pairs.coord
rc = params.Assemble_pairs_ref_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_ref_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_ref_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_ref_assemble_pairs.failed
fasta = params.Assemble_pairs_ref_assemble_pairs.fasta
nproc = params.Assemble_pairs_ref_assemble_pairs.nproc
alpha = params.Assemble_pairs_ref_assemble_pairs.alpha
maxerror = params.Assemble_pairs_ref_assemble_pairs.maxerror
minlen = params.Assemble_pairs_ref_assemble_pairs.minlen
maxlen = params.Assemble_pairs_ref_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_ref_assemble_pairs.scanrev
minident = params.Assemble_pairs_ref_assemble_pairs.minident
evalue = params.Assemble_pairs_ref_assemble_pairs.evalue
maxhits = params.Assemble_pairs_ref_assemble_pairs.maxhits
fill = params.Assemble_pairs_ref_assemble_pairs.fill
aligner = params.Assemble_pairs_ref_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_ref_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_ref_assemble_pairs.// dbexec
gap = params.Assemble_pairs_ref_assemble_pairs.gap
usearch_version = params.Assemble_pairs_ref_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_ref_assemble_pairs.assemble_reference
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi
	
	echo \$align_exec
	echo \$dbexec
	
	
	
	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc} >> out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g14_12_reads0_g2_0
 val mate from g_10_mate_g2_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g2_0_reads0_g_21
 set val(name), file("FS_*")  into g2_0_logFile1_g2_5
 set val(name), file("*_${method}-fail.fastq") optional true  into g2_0_reads22
 set val(name),file("out*") optional true  into g2_0_logFile3_g12_0

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed >> out_${R1}_FS.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

input:
 set val(name), file(log_file) from g2_0_logFile1_g2_5
 val mate from g_10_mate_g2_5

output:
 set val(name), file("*.tab")  into g2_5_logFile0_g2_7

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

input:
 val matee from g_10_mate_g2_7
 set val(name), file(log_files) from g2_5_logFile0_g2_7

output:
 file "*.rmd"  into g2_7_rMarkdown0_g2_9


shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_render_rmarkdown {

input:
 file rmk from g2_7_rMarkdown0_g2_9

output:
 file "*.html"  into g2_9_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_ref_parse_log_AP {

input:
 set val(name),file(log_file) from g14_12_logFile1_g14_15
 val mate from g_7_mate_g14_15

output:
 set val(name),file("*.tab")  into g14_15_logFile0_g14_19

script:
field_to_parse = params.Assemble_pairs_ref_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_ref_report_assemble_pairs {

input:
 set val(name),file(log_files) from g14_15_logFile0_g14_19
 val matee from g_7_mate_g14_19

output:
 file "*.rmd"  into g14_19_rMarkdown0_g14_22



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_ref_render_rmarkdown {

input:
 file rmk from g14_19_rMarkdown0_g14_22

output:
 file "*.html"  into g14_22_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality1_filter_seq_quality {

input:
 set val(name),file(reads) from g1_12_reads0_g19_0
 val mate from g_10_mate_g19_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g19_0_reads0_g_21
 set val(name), file("FS_*")  into g19_0_logFile1_g19_5
 set val(name), file("*_${method}-fail.fastq") optional true  into g19_0_reads22
 set val(name),file("out*") optional true  into g19_0_logFile33

script:
method = params.Filter_Sequence_Quality1_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality1_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality1_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality1_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality1_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed >> out_${R1}_FS.log
	"""
}


}


process macc_process_cat_align_and__ref {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /new.fasta$/) "final_reads/$filename"}
input:
 set val(name),  file(reads_align) from g19_0_reads0_g_21
 set val(name),  file(reads_ref) from g2_0_reads0_g_21

output:
 file "new.fasta"  into g_21_fasta_file00


script:
	
readArray_reads_align = reads_align.toString().split(' ')[0]	
readArray_reads_ref = reads_ref.toString().split(' ')[0]	

"""

sed -n '1~4s/^@/>/p;2~4p' ${readArray_reads_align} > output.fasta
sed -n '1~4s/^@/>/p;2~4p' ${readArray_reads_ref} > output1.fasta

sed -i 's/^>/>_ref/g' output1.fasta

cat  output1.fasta output.fastq > new.fasta
 
"""
}


process Filter_Sequence_Quality1_parse_log_FS {

input:
 set val(name), file(log_file) from g19_0_logFile1_g19_5
 val mate from g_10_mate_g19_5

output:
 set val(name), file("*.tab")  into g19_5_logFile0_g19_7, g19_5_logFile0_g12_0

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process make_report_pipeline_cat_all_file {

input:
 set val(name), file(log_file) from g1_12_logFile3_g12_0
 set val(name), file(log_file) from g2_0_logFile3_g12_0
 set val(name), file(log_file) from g19_5_logFile0_g12_0

output:
 set val(name), file("all_out_file.log")  into g12_0_logFile0_g12_2

script:
readArray = log_file.toString()

"""

echo $readArray
cat out* >> all_out_file.log
"""

}


process make_report_pipeline_report_pipeline {

input:
 set val(name), file(log_files) from g12_0_logFile0_g12_2

output:
 file "*.rmd"  into g12_2_rMarkdown0_g12_1


shell:

readArray = log_files.toString().split(' ')
R1 = readArray[0]

'''
#!/usr/bin/env perl


my $script = <<'EOF';


```{r, message=FALSE, echo=FALSE, results="hide"}
# Setup
library(prestor)
library(knitr)
library(captioner)

plot_titles <- c("Read 1", "Read 2")
if (!exists("tables")) { tables <- captioner(prefix="Table") }
if (!exists("figures")) { figures <- captioner(prefix="Figure") }
tables("count", 
       "The count of reads that passed and failed each processing step.")
figures("steps", 
        paste("The number of reads or read sets retained at each processing step. 
               Shown as raw counts (top) and percentages of input from the previous 
               step (bottom). Steps having more than one column display individual values for", 
              plot_titles[1], "(first column) and", plot_titles[2], "(second column)."))
```

```{r, echo=FALSE}
console_log <- loadConsoleLog(file.path(".","!{R1}"))
```

# Summary of Processing Steps

```{r, echo=FALSE}
count_df <- plotConsoleLog(console_log, sizing="figure")
```

`r figures("steps")`

```{r, echo=FALSE}
kable(count_df[c("step", "task", "total", "pass", "fail")],
      col.names=c("Step", "Task", "Input", "Passed", "Failed"),
      digits=3)
```

`r tables("count")`


EOF
	
open OUT, ">!{name}.rmd";
print OUT $script;
close OUT;

'''
}


process make_report_pipeline_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "outputparam/$filename"}
input:
 file rmk from g12_2_rMarkdown0_g12_1

output:
 file "*.html"  into g12_1_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality1_report_filter_Seq_Quality {

input:
 val matee from g_10_mate_g19_7
 set val(name), file(log_files) from g19_5_logFile0_g19_7

output:
 file "*.rmd"  into g19_7_rMarkdown0_g19_9


shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality1_render_rmarkdown {

input:
 file rmk from g19_7_rMarkdown0_g19_9

output:
 file "*.html"  into g19_9_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
