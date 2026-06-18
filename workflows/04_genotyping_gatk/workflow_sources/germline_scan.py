import sys
import vcf

def calculate_allele_balance(ad):
    if ad:
        ref = int(ad[0])
        alt = int(ad[1])
        total = ref + alt
        if total > 0:
            return alt / total
    return None

def main(vcf_file, father_id, mother_id, child_id, output_file):
    vcf_reader = vcf.Reader(filename=vcf_file)
    vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)

    for record in vcf_reader:
        # Check if the FILTER field is PASS or empty
        if record.FILTER is None or record.FILTER == [] or 'PASS' in record.FILTER:
            father_gt = record.genotype(father_id)['GT']
            mother_gt = record.genotype(mother_id)['GT']
            child_gt = record.genotype(child_id)['GT']
            child_ad = record.genotype(child_id)['AD']

            if father_gt == "0/0" and mother_gt == "0/0" and child_gt == "0/1":
#                print(child_ad)
#                allele_balance = calculate_allele_balance(child_ad)
#                if 0.3 <= allele_balance <= 0.7:
                    vcf_writer.write_record(record)

    vcf_writer.close()

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <vcf_file> <father_id> <mother_id> <child_id> <output_vcf_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

