
##add noise to data
if ! [ -x "$(command -v audio_degrader)" ]; then
  echo 'Error: audio_degrader is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v ffmpeg)" ]; then
  echo 'Error: ffmpeg is not installed.' >&2
  exit 1
fi

if [ $# -ne 3 ];then
        echo "./generate_noisy_dataset.sh <data_set_directory> <noise_file> <output_dir>"
        exit 1
fi

out_dir=$3

input_dir=$1
file_noise=$2

if [[ -f $input_dir ]]; then
    echo "$input_dir is a file. required a directory"
    exit 1
fi
if [[ -f $out_dir ]]; then
    echo "$out_dir is a file. required a directory"
    exit 1
fi
if [[ -d $file_noise ]]; then
    echo "$file_noise is a directory. required a file"
    exit 1
fi

if [ ! -d $out_dir ];then
	mkdir $out_dir
fi
for file in $input_dir/*.wav;do
	dir1="$(echo $file | cut -d'/' -f 1)"
	echo $dir1
	
	wavfilepath="$(basename $(dirname $file))/$(basename $file)"
	noisewavfilepath="$(basename $(dirname $file_noise))/$(basename $file_noise)"
	#echo $wavfilepath
	wavfile="$(echo ${wavfilepath##*/})"
	noisewavfile="$(echo ${noisewavfilepath##*/})"
	echo $wavfile
	echo $noisewavfile
	ffmpeg -i $file original_tmp.wav
	#length of the noise file
	time="$(soxi -D $noisewavfile)"
	#echo $time
	half="$(echo "$time/2" | bc -l)"
	#echo $half
	halff=`echo ${half%.*}`
	#echo $halff
	a=$((halff+1))
	#echo $a
	b=$((halff-1))
	#echo $b
	timef=`echo ${time%.*}`
	#randomly select a portion of noise
	start=`shuf -i 0-$a -n 1`
	#echo $start
	end=`shuf -i $b-$timef -n 1`
	#echo $end
	ffmpeg -i batman.wav -ss $start -to $end -c copy noise_tmp.wav
	audio_degrader -i original_tmp.wav -d mix,noise_tmp.wav//6 -o $out_dir/noise_$wavfile.wav
	#aplay degraded1.wav
	rm original_tmp.wav
	rm noise_tmp.wav
done
