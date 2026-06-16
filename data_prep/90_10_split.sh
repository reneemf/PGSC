#!/usr/bin/env bash

#SBATCH -J 90_10_split
#SBATCH --mem=25GB
#SBATCH --time=2:00:00
#SBATCH --partition=tier3q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/test_split/90_10_split_%A_%a.out
#SBATCH -e /home/<username>/slurm_outputs/test_split/90_10_split_%A_%a.err

now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"

# 90/10 split where white brit is train/test pop & white euro is validate
source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
pop_path="${JOINED_GENOS}${TRAIN_POP}/"
split_genos="${JOINED_GENOS}${TRAIN_POP}_split/"
pop_file="ukb_chr1-22_${TRAIN_POP}"

mkdir -p "$split_genos"

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

# randomly shuffle fam file
shuf --random-source=<(get_seeded_random 42) \
  "${pop_path}${pop_file}.fam" \
  -o "${pop_path}/${pop_file}_shuf.fam"

# split into 10 folds
split -d -n l/10 -a 1 \
  "${pop_path}${pop_file}_shuf.fam" \
  "${split_genos}${pop_file}_fold_"

# Generate testing and training sets
test_set="${split_genos}${pop_file}_test"
mv "${split_genos}${pop_file}_fold_0" "$test_set"
training_set="${split_genos}${pop_file}_train"

# Combine all training folds into one file
cat "${split_genos}${pop_file}_fold_"[1-9] > "$training_set"

# clean up intermediate folds
rm "${split_genos}${pop_file}_fold_"[1-9]
rm "${pop_path}${pop_file}_shuf.fam"

# sanity checks
echo "Original samples:  $(wc -l < "${pop_path}${pop_file}.fam")"
echo "Training samples:  $(wc -l < "$training_set")"
echo "Testing samples:   $(wc -l < "$test_set")"

# ensure no overlap
overlap=$(comm -12 \
  <(sort "$training_set") \
  <(sort "$test_set") | wc -l)

if [ "$overlap" -ne 0 ]; then
  echo "ERROR: Train/test overlap detected"
  exit 1
fi

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 

