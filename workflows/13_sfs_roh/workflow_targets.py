from workflow_templates import *
# Generating a list from "aa" to "cr" in Python
import string
import os

def generate_alphabet_list(start, end):
    alphabet = string.ascii_lowercase
    result = []

    # Finding the starting and ending positions
    start_index = alphabet.index(start[0]) * len(alphabet) + alphabet.index(start[1])
    end_index = alphabet.index(end[0]) * len(alphabet) + alphabet.index(end[1])

    # Generating the list
    for i in range(start_index, end_index + 1):
        first_char = alphabet[i // len(alphabet)]
        second_char = alphabet[i % len(alphabet)]
        result.append(first_char + second_char)

    return result


def generate_vcf_file(input_string):
    # Template for the file content
    template = """/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations_refined/genotype_trio/D_family1_1/D_family1_1_combine_{input}.vcf.gz
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations_refined/genotype_trio/D_family2_1/D_family2_1_combine_{input}.vcf.gz
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations_refined/genotype_trio/D_family3_1/D_family3_1_combine_{input}.vcf.gz
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations_refined/genotype_trio/D_family4_1/D_family4_1_combine_{input}.vcf.gz
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations_refined/genotype_trio/D_family5_1/D_family5_1_combine_{input}.vcf.gz
"""
    # Replace the placeholder with the input string
    file_content = template.format(input=input_string)
    # Path for the output file
    file_path = '/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/SFS/DUM/vcf_lists/DUM_{input}_vcfs.txt'.format(input = input_string)
    if os.path.exists(file_path):
        return f"File already exists: {file_path}"
    else:
        # Writing to the file
        with open(file_path, 'w') as file:
            file.write(file_content)
        return  f"File created: {file_path}"

def merge_and_filter(gwf,part):
    generate_vcf_file(part)
    vcfs = '/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/SFS/DUM/vcf_lists/DUM_{input}_vcfs.txt'.format(input = part)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/SFS/DUM/{input}".format(input = part)
    inds = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/SFS/DUM_inds.txt"
    outname = "DUM_{input}".format(input = part)
    gwf.target_from_template(
        name = "process_vcfs_{outname}".format(outname=outname),
        template = process_vcfs(path,vcfs,inds,outname)
    )
