#/usr/bin/python
__author__ = 'mcobb'

import sys
import os.path


def main():
	# Reads the file in from the command line. First file is the script, second is the vcf file,
	# and an option second is the summary file.
	num_of_files = len(sys.argv[1:])
	split_file = sys.argv[1]
	output = sys.argv[2]

	# If only a vcf file. Will be adjusted to automatically no through the command line input.
	if num_of_files == 2:
		# collect and import vcf file.
		# if os.path.isfile(split_file):
		# 	print "found file."
		# else:
		# 	print "vcf_path was not found"
		with open(split_file, 'r') as vcf_file:

			print "Results found at: ", output

			with open(output, 'wb') as new_file:
				for x in vcf_file:
					if x.startswith("Ld"):
						snp = x.split()
						# format = snp[8].split(':')
						# PL_index = [i for i, s in enumerate(format) if s.startswith('PL')]
						try:
							libs = snp[9:]
							for lib in libs:
								index = libs.index(lib)
								if lib == '.' or lib == './.:.:.:.:.':
									replace = '0/0:99:56:100,0:0,1992,1992'
									libs[index] = replace
							snp[9:] = libs
							line = '\t'.join(snp)
							new_file.write(line + "\n")
						except IndexError:
							# print snp
							pass

					else:
						new_file.write(x)
			new_file.close()
		vcf_file.close()

main()

		# path = os.path.dirname(os.path.abspath(__file__))
		# file_path = path + '/' + vcf_path.strip('.vcf') + '_replace.vcf'