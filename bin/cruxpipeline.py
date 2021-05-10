#!/usr/bin/env python3
# CREATE DATE: 9 May 2019
# AUTHOR: William Stafford Noble
import sys
import os
import subprocess

USAGE = """USAGE: crux-pipeline.py [options] <fasta> <ms2 directory>

  This program uses Crux to analyze a given set of MS/MS files.  The
  steps are
  
  1. Search using Tide or Comet.
  2. Process search results using Percolator.
  3. Aggregate counts for different isoforms.
  4. Compute NSAF values using spectral-counts.

  The required inputs are (1) a protein database in FASTA format and
  (2) a directory containing files of MS/MS spectra.  The latter can
  be in any format supported by ProteoWizard. On Windows, this can
  include RAW format; on Linux, RAW files must first be converted to
  some other format using msconvert.

  The output is a directory containing a summary file and a series of
  subdirectories.  The name of each subdirectory will be created based
  on the name of the input MS/MS file (minus its extension).  Each
  subdirectory will contain a series of intermediate files, as well as
  the spectral counts output, which is called spectral-counts.txt,
  containing protein-level NSAF scores.  In the end, all of the
  individual spectral-counts output files are summarized in a
  composite file called spectral-counts.txt.

  Note that with Tide, the program will also create a directory
  containing the database index. This will have the same name as the
  fasta file, minus the extension.  In subsequent runs of the
  pipeline, this index will be re-used if it already exists.

  This script assumes that the Crux binary is located in the path.

  Options:

   --search-engine tide-search|comet Specify the search engine 
                                     (default=tide-search)

"""

###############################################################################
PARAMETERS = """# comet_version 2016.01 rev. 2
# This file was created automatically by crux-pipeline.py.
# Everything following the '#' symbol is treated as a comment.

# This file is based on the Comet parameter file used by Momo, but it
# has been edited to include corresponding Tide parameters. --WSN

#database_name=swiss.HUMAN.20141123.amyloid.v2.concatenated.fasta
# Crux does not like having the database name in the parameter file. --WSN
decoy_search=1 # This was originally 0. --WSN
# 0=no (default), 1=concatenated search, 2=separate search
concat=T
decoy-format=shuffle

num_threads=0                        # 0=poll CPU to set num threads; else specify num threads directly (max 64)
num-threads=1 # Use 1 thread to avoid threading bug. --WSN

#
# masses
#
peptide_mass_tolerance=10.00
precursor-window=10
peptide_mass_units=2                 # 0=amu, 1=mmu, 2=ppm
precursor-window-type=ppm
mass_type_parent=1                   # 0=average masses, 1=monoisotopic masses
mass_type_fragment=1                 # 0=average masses, 1=monoisotopic masses
precursor_tolerance_type=0           # 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances
isotope_error=0                      # 0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)

#
# search enzyme
#
search_enzyme_number=1               # choose from list at end of this params file
enzyme=trypsin
num_enzyme_termini=1                 # 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific
digestion=partial-digest
allowed_missed_cleavage=2            # maximum value is 5; for enzyme search
missed-cleavages=2

#
# Up to 9 variable modifications are supported
# format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required>
#     e.g. 79.966331 STY 0 3 -1 0 0
#
variable_mod01=0.0 X 0 3 -1 0 0
variable_mod02=0.0 X 0 3 -1 0 0
variable_mod03=0.0 X 0 3 -1 0 0
variable_mod04=0.0 X 0 3 -1 0 0
variable_mod05=0.0 X 0 3 -1 0 0
variable_mod06=0.0 X 0 3 -1 0 0
variable_mod07=0.0 X 0 3 -1 0 0
variable_mod08=0.0 X 0 3 -1 0 0
variable_mod09=0.0 X 0 3 -1 0 0
max_variable_mods_in_peptide=5
require_variable_mod=0

#
# fragment ions
#
# ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions=1
# high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions=0
#
fragment_bin_tol=0.02              # binning to use on fragment ions
mz-bin-width=0.02
fragment_bin_offset=0.0              # offset position to start the binning (0.0 to 1.0)
mz-bin-offset=0.0
theoretical_fragment_ions=0          # 0=use flanking peaks, 1=M peak only
use-flanking-peaks=T
use_A_ions=0
use_B_ions=1
use_C_ions=0
use_X_ions=0
use_Y_ions=1
use_Z_ions=0
use_NL_ions=0                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks
use-neutral-loss-peaks=F

#
# output
#
output_sqtfile=1                     # 0=no, 1=yes  write sqt file
sqt-output=True
output_txtfile=0                     # 0=no, 1=yes  write tab-delimited txt file
txt-output=T # I include this for debugging. --WSN
output_pepxmlfile=1                  # 0=no, 1=yes  write pep.xml file
pepxml-output=T
output_percolatorfile=1              # 0=no, 1=yes  write Percolator tab-delimited input file
pin-output=T
output_outfiles=0                    # 0=no, 1=yes  write .out files
print_expect_score=1                 # 0=no, 1=yes to replace Sp with expect in out & sqt
num_output_lines=5                   # num peptide results to show
top-match=5
show_fragment_ions=0                 # 0=no, 1=yes for out files only

sample_enzyme_number=1               # Sample enzyme which is possibly different than the one applied to the search.
                                       # Used to calculate NTT & NMC in pepXML output (default=1 for trypsin).

#
# mzXML parameters
#
scan_range=0 0                       # start and scan scan range to search; 0 as 1st entry ignores parameter
precursor_charge=0 0                 # precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
override_charge=0                    # 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online
ms_level=2                           # MS level to analyze, valid are levels 2 (default) or 3
activation_method=ALL                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD

#
# misc parameters
#
digest_mass_range=600.0 5000.0       # MH+ peptide mass range to analyze
min-mass=600
max-mass=5000
num_results=100                      # number of search hits to store internally
skip_researching=1                   # for '.out' file output only, 0=search everything again (default), 1=don't search if .out exists
max_fragment_charge=3                # set maximum fragment charge state to analyze (allowed max 5)
max_precursor_charge=5               # set maximum precursor charge state to analyze (allowed max 9)
max-precursor-charge=5
nucleotide_reading_frame=0           # 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
clip_nterm_methionine=0              # 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine
clip-nterm-methionine=F
spectrum_batch_size=2000                # max. # of spectra to search at a time; 0 to search the entire scan range in one loop
decoy_prefix=decoy_                  # decoy entries are denoted by this string which is pre-pended to each protein accession
output_suffix=                       # add a suffix to output base names i.e. suffix "-C" generates base-C.pep.xml from base.mzXML input
mass_offsets=                        # one or more mass offsets to search (values substracted from deconvoluted precursor mass)

#
# spectral processing
#
minimum_peaks=10                     # required minimum number of peaks in spectrum to search (default 10)
min-peaks=10
minimum_intensity=0                  # minimum intensity value to read in
remove_precursor_peak=1              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)
remove-precursor-peak=T # This was originally False.
remove_precursor_tolerance=1.5       # +- Da tolerance for precursor removal
remove-precursor-tolerance=1.5
clear_mz_range=0.0 0.0               # for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range

#
# additional modifications
#

add_Cterm_peptide=0.0
add_Nterm_peptide=0.0
add_Cterm_protein=0.0
add_Nterm_protein=0.0

add_G_glycine=0.0000                 # added to G - avg.  57.0513, mono.  57.02146
add_A_alanine=0.0000                 # added to A - avg.  71.0779, mono.  71.03711
add_S_serine=0.0000                  # added to S - avg.  87.0773, mono.  87.03203
add_P_proline=0.0000                 # added to P - avg.  97.1152, mono.  97.05276
add_V_valine=0.0000                  # added to V - avg.  99.1311, mono.  99.06841
add_T_threonine=0.0000               # added to T - avg. 101.1038, mono. 101.04768
add_C_cysteine=57.021464             # added to C - avg. 103.1429, mono. 103.00918
add_L_leucine=0.0000                 # added to L - avg. 113.1576, mono. 113.08406
add_I_isoleucine=0.0000              # added to I - avg. 113.1576, mono. 113.08406
add_N_asparagine=0.0000              # added to N - avg. 114.1026, mono. 114.04293
add_D_aspartic_acid=0.0000           # added to D - avg. 115.0874, mono. 115.02694
add_Q_glutamine=0.0000               # added to Q - avg. 128.1292, mono. 128.05858
add_K_lysine=0.0000                  # added to K - avg. 128.1723, mono. 128.09496
add_E_glutamic_acid=0.0000           # added to E - avg. 129.1140, mono. 129.04259
add_M_methionine=0.0000              # added to M - avg. 131.1961, mono. 131.04048
add_O_ornithine=0.0000               # added to O - avg. 132.1610, mono  132.08988
add_H_histidine=0.0000               # added to H - avg. 137.1393, mono. 137.05891
add_F_phenylalanine=0.0000           # added to F - avg. 147.1739, mono. 147.06841
add_U_selenocysteine=0.0000          # added to U - avg. 150.3079, mono. 150.95363
add_R_arginine=0.0000                # added to R - avg. 156.1857, mono. 156.10111
add_Y_tyrosine=0.0000                # added to Y - avg. 163.0633, mono. 163.06333
add_W_tryptophan=0.0000              # added to W - avg. 186.0793, mono. 186.07931
add_B_user_amino_acid=0.0000         # added to B - avg.   0.0000, mono.   0.00000
add_J_user_amino_acid=0.0000         # added to J - avg.   0.0000, mono.   0.00000
add_X_user_amino_acid=0.0000         # added to X - avg.   0.0000, mono.   0.00000
add_Z_user_amino_acid=0.0000         # added to Z - avg.   0.0000, mono.   0.00000

#
# COMET_ENZYME_INFO _must_ be at the end of this parameters file
#
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P


"""

###############################################################################
def run_crux (arguments):
  command = ["crux"]
  command.extend(arguments.split())
  try:
    subprocess.check_call(command)
  except:
    sys.stderr.write("Error running command: %s\n" % " ".join(command))
    sys.exit(1)

###############################################################################
def extract_fasta_descriptions(fasta_file_name):
  """Extract from a FASTA file a dictionary in which the key is the
 sequence ID and value is the description."""

  return_value = {}

  fasta_file = open(fasta_file_name, "r")
  for line in fasta_file:
    line = line.rstrip()
    if (len(line) == 0):
      continue
    if (line[0] == ">"):

      # Parse the line.
      words = line[1:].split()
      identifier = words[0]
      description = " ".join(words[1:])

      # Store the description.
      return_value[identifier] = description
        
  fasta_file.close()
  sys.stderr.write("Read %d entries from %s.\n"
                   % (len(return_value), fasta_file_name))
  return(return_value)

###############################################################################
def aggregate_counts_files(counts_file_names, descriptions, aggregated_file_name):
  "Aggregate spectral counts from many files into one file."

  # Read all the spectral counts into memory.
  file_names = [] # List of basenames from files.
  identifiers = {} # Key = protein identifier, value = True
  spectral_counts = {} # Key = (file_name, identifier), value = count
  for counts_file_name in counts_file_names:

    # Use the name of the directory in which this file resides.
    file_name = os.path.basename(os.path.dirname(counts_file_name))
    file_names.append(file_name)

    counts_file = open(counts_file_name, "r")
    counts_file.readline() # Skip the header.
    num_counts = 0
    for line in counts_file:
      num_counts += 1
      try:
        (identifier, count) = line.rstrip().split("\t")
      except:
        sys.stderr.write("Error parsing line %d of %s (line=%s).\n"
                         % (num_counts, counts_file_name, line.rstrip()))
        sys.exit(1)
      identifiers[identifier] = True
      spectral_counts[(file_name, identifier)] = count
    sys.stderr.write("Read %d counts from %s.\n" % (num_counts, counts_file_name))
  sys.stderr.write("Found %d proteins in %d files.\n" % (len(identifiers), len(file_names)))

  # Make a table of spectral counts.
  if (aggregated_file_name == "-"):
    aggregated_file = sys.stdout
  else:
    aggregated_file = open(aggregated_file_name, "w")
  aggregated_file.write("Protein\tDescription")
  for file_name in file_names:
    aggregated_file.write("\t%s" % file_name)
  aggregated_file.write("\n")
  for identifier in identifiers.keys():
    aggregated_file.write("%s\t%s" % (identifier, descriptions[identifier]))
    for file_name in file_names:
      try:
        aggregated_file.write("\t%s" % spectral_counts[(file_name, identifier)])
      except:
        aggregated_file.write("\t")
    aggregated_file.write("\n")
  aggregated_file.close()        

###############################################################################
def remove_isoform (protein_id):
  """IDs of the form "<x>|<y>-<n>|<z>" are replaced with "<x>|<y>|<z>." """
  
  try:
    (first, second, third) = protein_id.split("|")
    try:
      (name, isoform) = second.split("-")
      return("%s|%s|%s" % (first, name, third))
    except:
      return(protein_id)
  except:
    return(protein_id)
  

###############################################################################
def rename_isoforms(percolator_output, no_isoforms_output):
  """Given a tab-delimited file containing one PSM per line, remove
     isoform numbers from the protein IDs."""

  percolator_file = open(percolator_output, "r")
  if (no_isoforms_output == "-"):
    no_isoforms_file = sys.stdout
  else:
    no_isoforms_file = open(no_isoforms_output, "w")

  # Read the header and identify the protein ID column.
  header = percolator_file.readline()
  no_isoforms_file.write(header)
  words = header.rstrip().split("\t")
  try:
    id_index = words.index("protein id")
  except:
    sys.stderr.write("Cannot find protein ID column in %s.\n"
                     % percolator_output)
    sys.stderr.write("%s\n" % "\n".join(words))
    sys.exit(1)

  # Traverse the file.
  line_num = 0
  for line in percolator_file:
    line_num += 1
    words = line.rstrip().split("\t")

    # Apply the transformation.
    new_protein_ids = []
    for protein_id in words[id_index].split(","):
      new_protein_id = remove_isoform(protein_id)
      if (not new_protein_id in new_protein_ids):
        new_protein_ids.append(new_protein_id)
    words[id_index] = ",".join(new_protein_ids)

    # Write out the new entry.
    no_isoforms_file.write("%s\n" % "\t".join(words))
  sys.stderr.write("Found %d entries in %s.\n" %
                   (line_num, percolator_output))
    
        
        
###############################################################################
# MAIN
###############################################################################
def main():
  global USAGE

  search_engine = "tide-search"

  # Parse the command line.
  sys.argv = sys.argv[1:]
  while (len(sys.argv) > 2):
      next_arg = sys.argv[0]
      sys.argv = sys.argv[1:]
      if (next_arg == "--search-engine"):
          if (sys.argv[0] == "comet") or (sys.argv[0] == "tide-search"):
              search_engine = sys.argv[0]
          else:
              sys.stderr.write("Unknown search engine (%s).\n"
                               % sys.argv[0])
              sys.exit(1)
          sys.argv = sys.argv[1:]
      else:
          sys.stderr.write("Unknown option (%s).\n" % next_arg)
          sys.exit(1)
  if (len(sys.argv) != 2):
      print(USAGE)
      sys.exit(1)
  fasta_file_name = sys.argv[0]
  msms_dir_name = sys.argv[1]

  # Create the parameter file.
  parameter_file_name = os.path.join(msms_dir_name, "parameters.txt")
  try:
    parameter_file = open(parameter_file_name, "w")
  except:
    sys.stderr.write("Cannot find directory %s.\n" % msms_dir_name)
    sys.exit(1)
  parameter_file.write(PARAMETERS)
  parameter_file.close()

  # All Crux commands share these arguments.
  shared_arguments = "--parameter-file %s" % parameter_file_name

  # If necessary, create the Tide index.
  if (search_engine == "tide-search"):
      index = os.path.splitext(os.path.basename(fasta_file_name))[0]
      if not os.path.exists(index):
        arguments = "tide-index %s --peptide-list T --output-dir %s %s %s" \
                    % (shared_arguments, index, fasta_file_name, index)
        run_crux(arguments)
  else:
      index = fasta_file_name

  # Loop over all the files in the given directory.
  counts_file_names = []
  for msms_file_name in os.listdir(msms_dir_name):

    # Construct the output directory name.
    (basename, extension) = os.path.splitext(os.path.basename(msms_file_name))
    if (extension == ".gz"):
      basename = os.path.splitext(basename)[0]
    # N.B. Could extend this list to other legal file types.
    elif ( (extension != ".raw") and (extension != ".RAW")
           and (extension != ".ms2") ):
      continue
    root = os.path.join(msms_dir_name, basename)
    if os.path.exists(root):
      sys.stderr.write("Error: %s exists.\n" % root)
      sys.exit(1)
    else:
      try:
        os.mkdir(root)
      except:
        sys.stderr.write("Error creating directory %s.\n" % root)
        sys.exit(1)
    sys.stderr.write("Writing results to %s.\n" % root)

    # Run the Crux commands.
    arguments = "pipeline %s --output-dir %s --search-engine %s %s %s" \
                % (shared_arguments, root, search_engine,
                   os.path.join(msms_dir_name, msms_file_name), index)
    run_crux(arguments)

    # Rename isoforms.
    percolator_output = os.path.join(root, "percolator.target.psms.txt")
    no_isoforms_output = os.path.join(root,
                                      "percolator.target.no-isoforms.psms.txt")
    rename_isoforms(percolator_output, no_isoforms_output)

    # Do spectral counting.
    arguments = "spectral-counts %s --find-peptides F --protein-database %s --output-dir %s %s" \
                % (shared_arguments, fasta_file_name, root,
                   no_isoforms_output)
    run_crux(arguments)
    counts_file_names.append(os.path.join(root, "spectral-counts.txt"))

  # Extract description lines from the FASTA.
  descriptions = extract_fasta_descriptions(fasta_file_name)

  # Aggregate all the spectral counts into one file.
  aggregate_counts_files(counts_file_names, descriptions,
                         os.path.join(msms_dir_name, "spectral-counts.txt"))

if __name__ == "__main__":
  main()
