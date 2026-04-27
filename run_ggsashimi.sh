#!/usr/bin/env bash

set -e
set -o pipefail

############################################
# PARSE COMMAND LINE ARGUMENTS
############################################

#Initialize variables
sample_file=""
coord_file=""
min_coverage=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_file)
      sample_file="$2"
      shift 2
      ;;
    --sample_file=*)
      sample_file="${1#*=}"
      shift
      ;;
    --coord_file)
      coord_file="$2"
      shift 2
      ;;
    --coord_file=*)
      coord_file="${1#*=}"
      shift
      ;;
    --min_coverage)
      min_coverage="$2"
      shift 2
      ;;
    --min_coverage=*)
      min_coverage="${1#*=}"
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--sample_file <file>] [--coord_file <file>]"
      echo ""
      echo "Options:"
      echo "  --sample_file FILE              table of samples to plot that includes columns `sample_name`, `cavatica_project`, `cram_name`, `group`"
      echo "  --coord_file FILE               table of regions to plot with columns `region` and `name`"
      echo "  --min_coverage VALUE            Minimum number of reads supporting a junction to be drawn [default=10]"
      echo "  -h, --help                      Display usage information"
      exit 0
      ;;
    *)
      >&2 echo "Error: Invalid argument '$1'"
      exit 1
      ;;
  esac
done

############################################
# MOUNT CAVATICA PROJECT
############################################

# Pull specified cavatica project from first row of sample file
cavatica_project=$(
awk -F'\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="cavatica_project") {
            col=i
            break
        }
    }
    next
}
NR==2 {
    print $col
    exit
}
' "$sample_file"
)

echo "cavatica_project=$cavatica_project"

# Create cavatica dir if it does not already exist
mkdir -p cavatica

## Mount cavatica project
sbfs mount --profile default --project $cavatica_project cavatica

echo "Waiting for sbfs mount to become active..."
until mountpoint -q cavatica; do
  sleep 1
done
echo "Mount is active."

cavatica_dir=cavatica/projects/$cavatica_project

echo "Waiting for cram files to become visible..."
until ls $cavatica_dir/GTE* >/dev/null 2>&1; do  
  sleep 1
done

echo "cavatica project successfully mounted. Proceeding..."

############################################
# CREATE template bam mapping file
############################################

# define output name
bammap_out="tmp/bammap.tsv"

# make tmp dir if it doesn't exist
mkdir -p tmp

# generate bam map file from sample info
awk -F'\t' -v OFS='\t' '
NR==1 {
    for (i=1;i<=NF;i++) {
        gsub(/\r/, "", $i)              # strip CR if present
        if ($i ~ /^sample_name$/) s=i
        if ($i ~ /^group$/) g=i
    }
    if (!s || !g) {
        print "ERROR: Could not find sample_name or group column" > "/dev/stderr"
        exit 1
    }
    next
}
{
    printf "%s\t%s\t%s\n", $s, "bams/" $s "_placeholder.bam", $g
}
' "$sample_file" > "$bammap_out"

echo "✅ Generated: $bammap_out"

############################################
# WRITE SAMPLE PROCESSING AND PLOTTING SCRIPT
############################################

echo "generating plot sashimi script based on input samples..."

output_script="plot_ggsashimi.sh"

echo "#!/usr/bin/env bash" > "$output_script"
echo "" >> "$output_script"

echo "" >> "$output_script"

echo "cavatica_dir=\"$cavatica_dir\"" >> "$output_script"
echo 'REGION=$1' >> "$output_script"
echo 'out=$2' >> "$output_script"
echo 'min_coverage=$3' >> "$output_script"
echo "" >> "$output_script"
echo "" >> "$output_script"

# LOOP THROUGH TSV AND ADD SAMTOOLS COMMANDS
echo "echo 'generating subsetted bams...'" >> "$output_script"

awk -F'\t' '
NR==1 {
    for (i=1;i<=NF;i++) {
        if ($i=="sample_name") s=i
        if ($i=="cram_name")   c=i
    }
    next
}
{
    print "samtools view -T refs/GRCh38.primary_assembly.genome.fa -b $cavatica_dir/"$c" ${REGION} -o tmp/bams/"$s"_$out.bam"
    print "samtools index tmp/bams/"$s"_$out.bam"
    print ""
}

' "$sample_file" >> "$output_script"

echo 'sed "s|placeholder|$out|g" tmp/bammap.tsv > tmp/bammap_$out.tsv' >>"$output_script"

echo "" >> "$output_script"

echo "echo 'generating sashimi plot...'" >> "$output_script"

echo 'python3 scripts/ggsashimi.py -b tmp/bammap_$out.tsv \
    -c ${REGION} \
    -g refs/gencode.v39.primary_assembly.annotation.gtf \
    -M $min_coverage -C 3 -O 3 \
    --alpha 1 --shrink --fix-y-scale \
    --overlay 3 --aggr mean_j \
    --base-size=14 \
    -R 350 --height=1.4 --width=7 \
    --ann-height 3 \
    -P examples/palette.txt \
    -o "output/sashimi_$out.pdf"' >> "$output_script"

chmod +x "$output_script"

echo "✅ Generated: $output_script"

############################################
# LOOP THROUGH REGIONS AND RUN PLOT SCRIPT
############################################

# Create tmp/bams dir if it does not already exist
mkdir -p tmp/bams
mkdir -p output

# loop through coord file rows to run ggsashimi plotting script
while read -r region name || [[ -n "$region" ]]; do

  echo "Processing $name region..."
  
  region=${region//$'\r'/}  
  name=${name//$'\r'/}

  bash $output_script "$region" "$name" $min_coverage; 

done < <(tail -n +2 "$coord_file")

# unmount cavatica project
sbfs unmount cavatica

# rm tmp files, plot_ggsashimi.sh script
rm -R tmp/*
rm plot_ggsashimi.sh

echo "sashimi plots generated ✅ "