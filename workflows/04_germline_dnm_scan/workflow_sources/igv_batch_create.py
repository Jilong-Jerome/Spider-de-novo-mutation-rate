import sys
def generate_igv_script(input_file, snapshot_directory, sp, ref_sp, bam_files_directory):
    geno_path = "/Users/au688344/GenomeDK/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome"
    igv_script = "igv_{sp}_script.txt".format(sp=sp)
    with open(input_file, 'r') as infile, open(igv_script, 'w') as outfile:
        for line in infile:
            cols = line.strip().split("\t")
            if len(cols) < 7: continue  # Skip incomplete lines
            
            chrom, pos, child_id, father_id, mother_id = cols[0], cols[1], cols[4], cols[5], cols[6]
            start_pos = str(int(pos) - 150)
            end_pos = str(int(pos) + 150)
            outfile.write("unload all\n") 
            outfile.write("new\n")
            outfile.write(f"genome {geno_path}/{ref_sp}_ncbi_chromosome.fa\n")
            outfile.write(f"snapshotDirectory {snapshot_directory}\n")
        
            # Load BAM files for the trio
            outfile.write(f"load {bam_files_directory}/{sp}/{mother_id}_final.bam\n")
            outfile.write(f"load {bam_files_directory}/{sp}/{father_id}_final.bam\n")
            outfile.write(f"load {bam_files_directory}/{sp}/{child_id}_final.bam\n")
            
            # Go to the specific region and take a snapshot
            outfile.write(f"goto {chrom}:{start_pos}-{end_pos}\n")
            outfile.write("sort position\n")
            outfile.write("squish\n")
            outfile.write("colorBy READ_STRAND\n")
            snapshot_name = f"{child_id}_{chrom}_{start_pos}_{end_pos}.png"
            outfile.write(f"snapshot {snapshot_name}\n")
        
        outfile.write("exit\n")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python generate_igv_script.py <input_file> <snapshot_directory> <species> <bam_files_directory>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    snapshot_directory = sys.argv[2]
    sp = sys.argv[3]
    if sp == "SAR":
        ref_sp = "SARA"
    elif sp == "BIC":
        ref_sp = "BI"
    elif sp == "TEN":
        ref_sp = "TENT"
    else:
        ref_sp = sp
    bam_files_directory = sys.argv[4]
    
    generate_igv_script(input_file, snapshot_directory, sp, ref_sp, bam_files_directory)
    print(f"IGV script has been generated. Run it in IGV using the -b option or through the GUI.")

