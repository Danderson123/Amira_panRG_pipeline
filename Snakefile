import os
import glob
import pandas as pd
import pyfastaq
import random
import statistics
from tqdm import tqdm
import json

#configfile: 'E_coli_config.yaml'
random.seed(10)
# get the list of all files in the directory
train_files = glob.glob(os.path.join(config["training_assembly_directory"], "*"))
test_files = []
input_files = train_files + test_files

# specify the output directory
output_dir = config["output_directory"]

rule all:
    input:
        os.path.join(output_dir, "AMR_supplemented_panRG.k15.w5.panidx.zip"),
        os.path.join(output_dir, "AMR_gene_headers_unified.txt"),
        os.path.join(output_dir, "AMR_alleles_unified.fa")


rule create_poppunk_input:
    input:
        input_files
    output:
        os.path.join(output_dir, "poppunk_input.txt")
    resources:
	    mem_mb=lambda wildcards, attempt: 10000 * attempt, threads=1
    run:
        # create an empty list with one sample per row
        samples = []
        # iterate through the input samples
        for sample_path in input:
            # get the sample name from the file path
            sample_name = os.path.basename(sample_path).split(".fa")[0]
            # tab separate the sample name and assembly path
            samples.append(f"{sample_name}\t{sample_path}")
        # write out the list
        with open(output[0], "w") as o:
            o.write("\n".join(samples))

rule run_poppunk:
    input:
        rules.create_poppunk_input.output
    output:
        directory(os.path.join(output_dir, "poppunk_output"))
    params:
        poppunk_database=config["poppunk_db_path"],
        skip_poppunk=config["skip_poppunk"]
    threads: config["threads"]
    resources:
        mem_mb=lambda wildcards, attempt: 30000 * attempt
    conda:
        "env/poppunk_env.yaml"
    shell:
        """
        if [ {params.skip_poppunk} != "True" ] && [ {params.skip_poppunk} != True ]; then
            poppunk_assign --run-qc --serial --max-zero-dist 1 --serial --write-references --max-merge 0 --db {params.poppunk_database} --query {input} --output {output} --threads {threads}
        else
            mkdir -p {output}
        fi
        """

checkpoint subsample_by_poppunk_cluster:
    input:
        poppunk_output=rules.run_poppunk.output,
        poppunk_input=rules.create_poppunk_input.output
    output:
        directory(os.path.join(output_dir, "subsampled_assemblies"))
    params:
        skip_poppunk=config["skip_poppunk"]
    run:
        # make the output directory
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # skip subsampling if specified
        if params.skip_poppunk == "True" or params.skip_poppunk == True:
            for f in input_files:
                shell(f"cp {f} {output[0]}")
                if ".gz" in f and not os.path.exists(f.replace(".gz", "")):
                    to_decompress = os.path.join(output[0], os.path.basename(f))
                    shell(f"gunzip {to_decompress}")
                if ".fna" in f:
                    name = os.path.join(output[0], os.path.basename(f).replace(".gz", "").replace(".fna", ".fa"))
                    shell(f"mv {os.path.join(output[0], os.path.basename(f)).replace('.gz', '')} {name}")
        else:
            # import the poppunk clusters
            cluster_assignments = pd.read_csv(os.path.join(input[0], "poppunk_output_clusters.csv"))
            # convert the poppunk input into a dictionary
            with open(input[1]) as i:
                content = i.read().split("\n")
            sample_paths = {}
            for row in content:
                tab_split = row.split("\t")
                sample_paths[tab_split[0].replace(".", "_")] = tab_split[1]
            # select 1 sample per poppunk cluster
            cluster_assignments = cluster_assignments[cluster_assignments['Taxon'].isin(sample_paths)]
            # get a list of the unique clusters
            clusters = list(set(cluster_assignments["Cluster"]))
            seen_samples = set()
            for c in tqdm(clusters):
                cluster_df = cluster_assignments[cluster_assignments['Cluster'] == c]
                assert not len(cluster_df) == 0
                sample_index = random.choice(cluster_df.index)  # Use random.choice for index
                sample_name = cluster_assignments.loc[sample_index, "Taxon"].replace("_query", "")
                assert not sample_name in seen_samples
                sample_path = sample_paths[sample_name]
                new_file_name = os.path.join(output[0], os.path.basename(sample_path))
                shell(f"cp {sample_path} {output[0]}")
                if ".gz" in new_file_name and not os.path.exists(new_file_name.replace(".gz", "")):
                    shell(f"gunzip {new_file_name} && rm -rf {new_file_name}")
                if "#" in new_file_name:
                    shell(f"mv {new_file_name.replace('.gz', '')} {new_file_name.replace('.gz', '').replace('#', '')}")
                if ".fna" in new_file_name:
                    shell(f"mv {new_file_name.replace('.gz', '')} {new_file_name.replace('.gz', '').replace('.fna', '.fa')}")
                seen_samples.add(sample_name)

rule make_AMR_gff:
    output:
        AMR_alleles=os.path.join(output_dir, "AMR_alleles.gff")
    params:
        ariba_path=config["ariba_path"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 10000, threads=1
    run:
        # get the output directory name
        output_dir = os.path.dirname(output[0])
        if not (os.path.exists(os.path.join(output_dir, "AMR_alleles.fa")) and os.path.exists(os.path.join(output_dir, "AMR_headers.txt"))):
            # make the output directory
            if not os.path.exists(os.path.join(output_dir, "ariba_downloads")):
                os.mkdir(os.path.join(output_dir, "ariba_downloads"))
            # download all AMR reference alleles using Ariba
            fasta_files = []
            for database in tqdm(["ncbi"]):#["argannot",
                            # "card",
                            # "resfinder",
                            # "srst2_argannot",
                            # "vfdb_core",
                            # "vfdb_full",
                            # "virulencefinder",
                            # "ncbi"]):
                db_output = os.path.join(output_dir, "ariba_downloads", database)
                #shell(f"singularity run {params.ariba_path} getref {database} {db_output}")
                shell(f"cp /hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG/ariba_downloads/* {os.path.dirname(db_output)}")
                fasta_files.append(db_output + ".fa")
            # cat all the reference AMR alleles
            all_alleles = []
            all_headers = []
            # iterate through the fasta files
            seen_headers = {}
            amr_calls = {}
            for f in tqdm(fasta_files):
                with open(f.replace(".fa", ".tsv")) as i:
                    rows = i.read().split("\n")
                for row in rows:
                    if row != " " and row != "":
                        amr_calls[row.split("\t")[0]] = row.split("\t")[5]
                # import the fasta file
                with open(f) as i:
                    content = i.read().split(">")[1:]
                # iterate through each entry in the fasta
                for sequence in content:
                    if not (sequence == "" or sequence == " "):
                        linesplit = sequence.split("\n")
                        header = linesplit[0]
                        if not header in seen_headers:
                            seen_headers[header] = 1
                        else:
                            seen_headers[header] += 1
                            header = header + "_" + str(seen_headers[header])
                        seq = "".join(linesplit[1:]).upper()
                        all_alleles.append(">" + header + "\n" + seq)
                        all_headers.append(header)
            # write out the files
            with open(os.path.join(output_dir, "AMR_alleles.fa"), "w") as o:
                o.write("\n".join(all_alleles))
            with open(os.path.join(output_dir, "AMR_headers.txt"), "w") as o:
                o.write("\n".join(all_headers))
            with open(os.path.join(output_dir, "AMR_calls.json"), "w") as o:
                o.write(json.dumps(amr_calls))
        # convert the FASTA to a GFF
        shell(f'python3 /nfs/research/zi/dander/general_bioinformatics_scripts/convert_fasta_to_gff.py {os.path.join(output_dir, "AMR_alleles.fa")} {output[0]}')

rule make_plasmid_gene_gff:
    input:
        config["plasmid_genes_path"]
    output:
        plasmid_alleles=os.path.join(output_dir, "plasmid_alleles.gff")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 10000, threads=1
    shell:
       "python3 /nfs/research/zi/dander/general_bioinformatics_scripts/convert_fasta_to_gff.py {input} {output}"

rule run_bakta:
    input:
        os.path.join(output_dir, "subsampled_assemblies", "{sample}.fa")
    output:
        directory(os.path.join(output_dir, "bakta_annotated_assemblies", '{sample}'))
    params:
        bakta_db=config["bakta_db_path"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 20000 * attempt, threads=1
    shell:
        "singularity run software/bakta.sif --db {params.bakta_db} --threads {threads} --output {output} {input}"
    #run:
    #   to_copy = os.path.join("/hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG_l_0/bakta_assemblies",
    #                   os.path.basename(output[0]))
    #   out = os.path.dirname(output[0])
    #   shell(f"cp -r {to_copy} {out}")

def get_samples():
    # This function should return a list of all samples that are processed by 'run_bakta'
    return [os.path.basename(f).replace(".fa", "") for f in glob.glob(os.path.join(checkpoint_output, "*.fa*"))]

def get_bakta_files(wildcards):
    checkpoint_output = checkpoints.subsample_by_poppunk_cluster.get(**wildcards).output[0]
    return expand(os.path.join(output_dir, "bakta_annotated_assemblies", "{sample}"), sample=[os.path.basename(f).replace(".fa", "") for f in glob.glob(os.path.join(checkpoint_output, "*.fa*"))])

rule remove_test_annotations:
    input:
        get_bakta_files
    output:
        touch(os.path.join(output_dir, "annotations_removed.done"))
    params:
        test_dir=config["testing_assembly_directory"]
    run:
        # list all gff files output by bakta
        gff_files = [f for f in glob.glob(os.path.join(os.path.dirname(input[0]), "*", "*.gff*")) + [input[-1]] if "edited" not in f]
        # list the names of the test samples
        test_samples = [os.path.basename(f).replace(".fa", "") for f in glob.glob(os.path.join(params.test_dir, "*"))]
        for g in tqdm(gff_files):
            for t in test_samples:
                if t in g:
                    print(g)
                    with open(g) as i:
                        features, sequence = i.read().split("##FASTA")
                    new_features = []
                    for l in features.split("\n"):
                        if l == "" or l == "\n":
                            continue
                        if l.startswith("#") or "CDS" not in l:
                            new_features.append(l)
                        else:
                            contig, source, cat, start, end, dot1, strand, dot2, annotations = l.split("\t")
                            split_annotations = annotations.split(";")
                            new_split_annotations = []
                            for i in split_annotations:
                                if "Name=" in i or "gene=" in i:
                                    continue
                                new_split_annotations.append(i)
                            new_features.append(f"{contig}\t{source}\t{cat}\t{start}\t{end}\t{dot1}\t{strand}\t{dot2}\t{';'.join(new_split_annotations)}")
                    with open(g, "w") as o:
                        o.write("\n".join(new_features + ["##FASTA", sequence]))

rule remove_short_annotations:
    input:
        get_bakta_files,
        rules.remove_test_annotations.output,
    output:
        touch(os.path.join(output_dir, "short_alleles_removed.done"))
    run:
        # list all gff files output by bakta
        gff_files = [f for f in glob.glob(os.path.join(os.path.dirname(input[0]), "*", "*.gff*"))]
        for g in tqdm(gff_files):
            with open(g) as i:
                features, sequence = i.read().split("##FASTA")
            new_features = []
            for l in features.split("\n"):
                if l == "" or l == "\n":
                    continue
                if l.startswith("#") or "CDS" not in l:
                    new_features.append(l)
                else:
                    contig, source, cat, start, end, dot1, strand, dot2, annotations = l.split("\t")
                    if int(end) - int(start) > 249:
                        new_features.append(f"{contig}\t{source}\t{cat}\t{start}\t{end}\t{dot1}\t{strand}\t{dot2}\t{annotations}")
            new_sequence = []
            for l in sequence.split("\n"):
                if l != "" and l != " ":
                    new_sequence.append(l)
            with open(g, "w") as o:
                o.write("\n".join(new_features + ["##FASTA", "\n".join(new_sequence)]))

rule list_bakta_gffs:
    input:
        get_bakta_files,
        rules.make_AMR_gff.output,
        rules.remove_test_annotations.output,
        rules.remove_short_annotations.output,
        rules.make_plasmid_gene_gff.output
    output:
        os.path.join(output_dir, "panaroo_input.txt")
    run:
        # list all gff files output by bakta
        gff_files = glob.glob(os.path.join(os.path.dirname(input[0]), "*", "*.gff*")) + [input[-4]] + [input[-1]]
        with open(output[0], 'w') as file_out:
            for gff in gff_files:
                file_out.write(gff + "\n")

rule run_panaroo:
    input:
        rules.list_bakta_gffs.output
    output:
        directory(os.path.join(output_dir, "panaroo_output"))
    params:
        identity=config["panaroo"]["identity"],
        len_dif_percent=config["panaroo"]["len_dif_percent"],
        length_outlier_support_proportion=config["panaroo"]["length_outlier_support_proportion"]
    threads: config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 40000 * attempt, threads=config["threads"]
    shell:
        "panaroo --clean-mode sensitive --refind-mode off --remove-invalid-genes -c {params.identity} --len_dif_percent {params.len_dif_percent} --length_outlier_support_proportion {params.length_outlier_support_proportion} --merge_paralogs -i {input} -o {output} --threads {threads}"

rule get_panaroo_alignments:
    input:
        rules.run_panaroo.output
    output:
        touch(os.path.join(output_dir, "panaroo_alignments.done"))
    threads: config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 40000 * attempt, threads=config["threads"]
    shell:
        #"panaroo-msa -o {input} -a pan --aligner mafft -t {threads}"
        #"python3 /hps/nobackup/iqbal/dander/amira_panRG_pipeline/software/panaroo/panaroo-msa-runner.py -o {input} -a pan --aligner mafft -t {threads}"
        "mv /hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG_plasmid_and_PAD_genes/panaroo_output/previous_aligned_gene_sequences /hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG_plasmid_and_PAD_genes/panaroo_output/aligned_gene_sequences"
        #"touch {output}"

checkpoint qc_panaroo_alignments:
    input:
        rules.get_panaroo_alignments.output
    output:
        directory(os.path.join(output_dir, "qced_unaligned_gene_sequences"))
    threads: 1
    resources:
	    mem_mb=lambda wildcards, attempt: 40000, threads=1
    params:
        test_dir=config["testing_assembly_directory"]
    run:
        def clean_gene(gene):
            chars_to_remove = set([".", "|", "(", ")", "-", "*", "+", "#", ":", "=", "/", ",","'"])
            cleaned_gene = "".join(char for char in gene if char not in chars_to_remove)
            return cleaned_gene

        # Ensure the output directory exists
        os.makedirs(output[0], exist_ok=True)
        input_dir = os.path.join(os.path.dirname(input[0]), "panaroo_output", "aligned_gene_sequences")
        # load the panaroo gene data csv
        gene_data = pd.read_csv(os.path.join(os.path.dirname(input_dir), "gene_data.csv"))
        # load the gene presence absence csv
        gene_presence_absence = pd.read_csv(os.path.join(os.path.dirname(input_dir), "gene_presence_absence.csv"))
        # get a list of locus IDs that correspond to transposases
        transposase_gene_annotation_ids = set()
        samples = list(gene_presence_absence.columns.values)[3:-1]
        for index, row in gene_presence_absence.iterrows():
            if isinstance(row["Annotation"], float):
                continue
            if "transposase" in row["Annotation"] or "Transposase" in row["Annotation"] or ("IS" in row["Annotation"] and "element" in row["Annotation"]):
                if isinstance(row["AMR_alleles"], float):
                    for samp in samples:
                        if not isinstance(row[samp], float):
                            if ";" in row[samp]:
                                for annotation_id in row[samp].split(";"):
                                    transposase_gene_annotation_ids.add(annotation_id.replace("_pseudo", ""))
                            else:
                                transposase_gene_annotation_ids.add(row[samp].replace("_pseudo", ""))
        # collect the geneIds of the AMR alleles
        gene_ids = {}
        transposase_genes = set()
        for index, row in gene_data.iterrows():
            if row["gff_file"] == "AMR_alleles":
                gene_ids[row["clustering_id"]] = row["scaffold_name"]
            assert not ";" in row["annotation_id"]
            if row["annotation_id"] in transposase_gene_annotation_ids:
                transposase_genes.add(row["clustering_id"])
            try:
                if "transposase" in row["description"] or "Transposase" in row["description"] or ("IS" in row["description"] and "element" in row["description"]):
                    transposase_genes.add(row["clustering_id"])
            except:
                pass
        # List and process alignments
        junk_alleles = []
        alleles = {}
        lengths = {}
        names = {}
        for a in tqdm(glob.glob(os.path.join(input_dir, "*"))):
            reader = pyfastaq.sequences.file_reader(a)
            alleles[a], lengths[a] = [], []

            # Process each sequence
            for sequence in reader:
                original_seq = str(sequence.seq)
                sequence.seq = str(sequence.seq).replace("-", "").upper()
                sample_id, cluster_id = sequence.id.split(";")
                if cluster_id in transposase_genes:
                    continue
                if len(sequence.seq) % 3 == 0 and "*" not in sequence.translate()[:-1]:
                    if not len(sequence.seq) > 249:
                        continue
                    alleles[a].append(f">{sequence.id}\n{sequence.seq}")
                    lengths[a].append(len(sequence.seq))
                    if cluster_id in gene_ids:
                        if not a in names:
                            names[a] = set()
                        names[a].add(clean_gene(gene_ids[cluster_id]))
                else:
                    junk_alleles.append(f">{sequence.id}\n{original_seq}")
        amr_allele_file_mapping = {}
        for a in tqdm(glob.glob(os.path.join(input_dir, "*"))):
            # Filter alleles based on length criteria
            if not lengths[a] == []:
                median_length = statistics.median(lengths[a])
                #length_tolerance = (0.8 * median_length, 1.2 * median_length)
                length_tolerance = (0, 5 * median_length)
                filtered_alleles = [allele for allele, length in zip(alleles[a], lengths[a]) if length_tolerance[0] <= length <= length_tolerance[1]]
                if not filtered_alleles == []:
                    if not all("refound" in allele for allele in filtered_alleles):
                        # Write to new file if any allele is ≥250bp
                        new_file_name = os.path.join(output[0], os.path.basename(a).replace("~~~", ".")).replace(".aln.fas", ".fasta")
                        if a in names:
                            amr_allele_file_mapping[new_file_name] = list(names[a])
                        with open(new_file_name, "w") as o:
                            if len(filtered_alleles) == 1:
                                allele = filtered_alleles[0].split("\n")
                                allele[1] = allele[1].replace("-", "")
                                o.write("\n".join(allele))
                            else:
                                o.write("\n".join(filtered_alleles))
        test_annotations = []
        # list the names of the test samples
        test_samples = set([os.path.basename(f).replace(".fa", "") for f in glob.glob(os.path.join(params.test_dir, "*"))])
        # iterate through the gced sequences
        for a in tqdm(glob.glob(os.path.join(output[0], "*"))):
            gene_name = os.path.basename(a).replace(".fasta", "")
            reader = pyfastaq.sequences.file_reader(a)
            # check if this fasta contains AMR alleles
            contains_AMR = False
            if a in amr_allele_file_mapping:
                contains_AMR = True
            # Process each sequence
            modified_sequences = []
            for sequence in reader:
                sample_id, cluster_id = sequence.id.split(";")
                if contains_AMR is False:
                    if sample_id in test_samples:
                        test_annotations.append(f">{gene_name};{sample_id};{cluster_id}\n{str(sequence.seq)}")
                    else:
                        modified_sequences.append(f">{sample_id};{cluster_id}\n{str(sequence.seq)}")
                if contains_AMR is True:
                    if "AMR_alleles" in sample_id:
                        modified_sequences.append(f">{sample_id};{cluster_id}\n{str(sequence.seq)}")
            if len(modified_sequences) != 0:
                with open(a, "w") as o:
                    o.write("\n".join(modified_sequences))
            else:
                os.remove(a)
        with open(os.path.join(os.path.dirname(input_dir), "failed_qc.fasta"), "w") as o:
            o.write("\n".join(junk_alleles))
        with open(os.path.join(os.path.dirname(input_dir), "test_sample_gene_annotations.fasta"), "w") as o:
            o.write("\n".join(test_annotations))
        import json
        with open(os.path.join(os.path.dirname(input_dir), "AMR_allele_to_Panaroo_COG_mapping.json"), "w") as o:
            o.write(json.dumps(amr_allele_file_mapping))

rule align_qced_sequences:
    input:
        unaligned=os.path.join(output_dir, "qced_unaligned_gene_sequences", "{gene}.fasta")
    output:
        aligned=os.path.join(output_dir, "qced_aligned_gene_sequences", "{gene}.fasta")
    threads: 2  # Adjust based on your system's capabilities
    resources:
	    mem_mb=lambda wildcards, attempt: 10000 * attempt, threads=2
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        with open(input[0]) as i:
            content = i.read().split(">")[1:]
        if len(content) > 1:
            shell("mafft --auto --thread {threads} {input[0]} > {output[0]}")
        else:
            shell("cp {input[0]} {output_dir}")

def get_processed_files(wildcards):
    checkpoint_output = checkpoints.qc_panaroo_alignments.get(**wildcards).output[0]
    return expand(os.path.join(output_dir, "qced_aligned_gene_sequences", "{gene}.fasta"), 
        gene=[os.path.splitext(os.path.basename(f))[0] 
        for f in glob.glob(os.path.join(output_dir, "qced_unaligned_gene_sequences", "*.fasta"))])

rule make_prg:
    input:
        get_processed_files
    output:
        os.path.join(output_dir, "AMR_supplemented_panRG.prg.fa")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 40000 * attempt, threads=16
    run:
        input_dir = os.path.dirname(input[0])
        prefix = output[0].replace(".prg.fa", "")
        shell("make_prg from_msa -F -i {input_dir} --output-prefix {prefix} --threads 32")

rule build_pandora_index:
    input:
        prg=rules.make_prg.output
    output:
        os.path.join(output_dir, "AMR_supplemented_panRG.k15.w5.panidx.zip")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 30000 * attempt, threads=16
    params:
        pandora="/hps/nobackup/iqbal/dander/Escherichia_coli_panRG_c_0.8_l_0_train_AMR_alleles_removed/software/pandora-linux-precompiled-v0.12.0-alpha.0",
        kmer_size=15,
        window_size=5
    shell:
        "{params.pandora} index -w {params.window_size} -k {params.kmer_size} -t 32 -o {output} {input.prg}"

rule make_processed_amr_headers:
    input:
        alignments=get_processed_files,
        amr_headers=rules.make_AMR_gff.output
    output:
        os.path.join(output_dir, "AMR_gene_headers_unified.txt")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 10000, threads=1
    run:
        def clean_gene_list(gene_list):
            chars_to_remove = set([".", "|", "(", ")", "-", "*", "+", "#", ":", "=", "/", ",","'"])
            cleaned_genesOfInterest = []
            for g in gene_list:
                cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
                cleaned_genesOfInterest.append(cleaned_gene)
            return set(cleaned_genesOfInterest)

        import json
        # import the json saying where the amr alleles are
        with open(os.path.join(os.path.dirname(input.amr_headers[0]), "panaroo_output/AMR_allele_to_Panaroo_COG_mapping.json")) as i:
            amr_mapping = json.load(i)
        cleaned = {}
        for key in amr_mapping:
            cleaned[os.path.basename(key)] = amr_mapping[key]
        amr_mapping = cleaned
        amr_genes = set()
        with open(os.path.join(os.path.dirname(input.amr_headers[0]), "AMR_headers.txt")) as i:
            genes_of_interest = clean_gene_list(i.read().split("\n"))
        for a in tqdm(input.alignments):
            if os.path.basename(a) in amr_mapping:
                gene_name = os.path.basename(a).replace(".fasta", "")
                amr_genes.add(gene_name)
        with open(output[0], "w") as o:
            o.write("\n".join(list(amr_genes)))

rule make_processed_amr_fasta:
    input:
        alignments=get_processed_files,
        amr_headers=rules.make_AMR_gff.output
    output:
        os.path.join(output_dir, "AMR_alleles_unified.fa")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 10000, threads=1
    run:
        def clean_gene_list(gene_list):
            chars_to_remove = set([".", "|", "(", ")", "-", "*", "+", "#", ":", "=", "/", ",","'"])
            cleaned_genesOfInterest = []
            for g in gene_list:
                cleaned_gene = "".join(char for char in g if char not in chars_to_remove)
                cleaned_genesOfInterest.append(cleaned_gene)
            return set(cleaned_genesOfInterest)

        import json
        # import the json saying where the amr alleles are
        with open(os.path.join(os.path.dirname(input.amr_headers[0]), "panaroo_output/AMR_allele_to_Panaroo_COG_mapping.json")) as i:
            amr_mapping = json.load(i)
        # import the reference alleles
        reader = pyfastaq.sequences.file_reader(os.path.join(os.path.basename(output_dir), "AMR_alleles.fa"))
        # Process each sequence
        ref_alleles = {}
        for sequence in reader:
            cleaned_name = list(clean_gene_list([sequence.id]))[0]
            assert cleaned_name not in ref_alleles
            ref_alleles[cleaned_name] = {"identifier": sequence.id, "sequence": str(sequence.seq)}
        # collect the non-ref alleles
        unified_fasta_content = []
        for fasta_file in amr_mapping:
            # get the cluster name
            unified_name = os.path.basename(fasta_file).replace(".fasta", "")
            for ref_allele_name in amr_mapping[fasta_file]:
                unified_fasta_content.append(f">{unified_name};{ref_alleles[ref_allele_name]['identifier']}\n{ref_alleles[ref_allele_name]['sequence']}")
            # load the fasta file
            reader = pyfastaq.sequences.file_reader(fasta_file)
            allele_count = 1
            #for sequence in reader:
            #    if "AMR_alleles" not in sequence.id:
            #        unified_fasta_content.append(f">{unified_name};SUPP_{unified_name}_{allele_count}\n{str(sequence.seq)}")
            #        allele_count += 1
        with open(output[0], "w") as o:
            o.write("\n".join(unified_fasta_content))

