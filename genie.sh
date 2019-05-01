ORIG_DATA="data/$2.fa"
PP_DATA="data/pp_$2.fa"
BWT_FILE="data/$2.bwt"
LUT_FILE="data/$2.lut"
SMALL=21
BIG=151
SMALL_QUERY="data/$2_small_1k.fa"
BIG_QUERY="data/$2_big_1k.fa"
LUT_OUTPUT_FILE="data/$2_output.lut"
RMI_OUTPUT_FILE="data/$2_output.rmi"
OUTPUT_FILE="output_$2"

if [ "$1" = "clean" ]; then
  echo "cleaning..."
  rm $PP_DATA
  rm $BWT_FILE
  rm $SMALL_QUERY
  rm $BIG_QUERY
  rm build_bwt_index/bwtbuildbinarybwt
  rm build_bwt_index/bwt-build-with-cp.8bit.32
  rm naive/bwt-naive
  rm opt/bwt-match-checkpoint.8bit.32
  rm $OUTPUT_FILE
  echo "...done cleaning"
fi

if [ "$1" == "run_rmi" ] || [ "$1" = "run_lut" ] || [ "$1" = "run_opt" ] ||  [ "$1" = "run_naive" ] ;  then
  echo "preprocessing..."
  cp $ORIG_DATA $PP_DATA
  echo "...done copying"
  gsed -i 's/N//g' $PP_DATA
  echo "...done 1"
  gsed -i '/^$/d'  $PP_DATA
  echo "...done 2"
  gsed -i 's/a/A/g' $PP_DATA
  echo "...done a"
  gsed -i 's/c/C/g' $PP_DATA
  echo "...done c"
  gsed -i 's/g/G/g' $PP_DATA
  echo "...done g"
  gsed -i 's/t/T/g' $PP_DATA
  echo "...done t"
  gsed -i '/>/d'  $PP_DATA
  gsed -i '/^$/d'  $PP_DATA
  gsed -i 's/n//g' $PP_DATA
  gsed -i '/^$/d' $PP_DATA
  gsed -i '1s/^/>refseq\n/' $PP_DATA
  echo "...done preprocessing"
fi

if [ "$1" == "run_rmi" ] || [ "$1" = "run_lut" ] || [ "$1" = "run_opt" ] ||  [ "$1" = "run_naive" ] ;  then
  echo "compiling..."
  g++ -std=c++11 -I seqan-library-2.1.1/include/ opt/run_lut.c build_bwt_index/file.cpp -o opt/run_lut
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with run_lut"
  g++ -std=c++11 -I seqan-library-2.1.1/include/ build_bwt_index/bwtbuildbinarybwt.cpp build_bwt_index/file.cpp -o build_bwt_index/bwtbuildbinarybwt
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with bwtbuildbinarybwt"
  g++ -std=c++11 -I. -I seqan-library-2.1.1/include/ build_bwt_index/bwt-build-with-cp.cpp build_bwt_index/file.cpp -o build_bwt_index/bwt-build-with-cp.8bit.32
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with bwt-build-with-cp"
  clang++ -Xpreprocessor -fopenmp -lomp -O1 -DPRINT_OUTPUT naive/bwt-naive.c build_bwt_index/file.cpp -o naive/bwt-naive
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with bwt-naive"
  clang++ -Xpreprocessor -fopenmp -lomp -O3 -v -march=haswell opt/bwt-match-checkpoint.c build_bwt_index/file.cpp -I. -o opt/bwt-match-checkpoint.8bit.32
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with bwt-match-checkpoint"
  clang++ -Xpreprocessor -fopenmp -lomp -O3 -v -march=haswell opt/bwt-clean-match-checkpoint.c build_bwt_index/file.cpp -I. -o opt/bwt-clean-match-checkpoint.8bit.32
  if [ $? -ne 0 ]; then
    exit 0
  fi
  echo ">>done with bwt-clean-match-checkpoint"

  echo "...done compiling"
fi

if [ "$1" == "run_rmi" ] || [ "$1" = "run_lut" ] || [ "$1" = "run_opt" ] || [ "$1" = "run_naive" ];  then
  echo "building compressed FM Index..."
  time ./build_bwt_index/bwtbuildbinarybwt $PP_DATA 1 $BWT_FILE $LUT_FILE
  time ./build_bwt_index/bwt-build-with-cp.8bit.32 $BWT_FILE $PP_DATA
  echo "...done building compressed FM Index..."
fi

if [ "$1" == "run_rmi" ] || [ "$1" = "run_lut" ] || [ "$1" = "run_opt" ] ||  [ "$1" = "run_naive" ] ;  then
  echo "creating queries..."
  python2.7 ./data/selectQueries.anyLength.py $PP_DATA 1000 $SMALL > $SMALL_QUERY
  python2.7 ./data/selectQueries.anyLength.py $PP_DATA 1000 $BIG > $BIG_QUERY
  echo "...done creating queries"
fi

if [ "$1" = "run_lut" ];  then
  echo "running lut .... "
  time ./opt/run_lut $LUT_FILE $BWT_FILE $PP_DATA $SMALL_QUERY $SMALL $LUT_OUTPUT_FILE
  echo " .... done running lut"
fi

if [ "$1" = "run_rmi" ]; then
  echo "running rmi.... "
  time python3 ./learned_index/train.py $BWT_FILE $PP_DATA $SMALL_QUERY $SMALL $RMI_OUTPUT_FILE
  echo " .... done running rmi"
fi

if [ "$1" = "run_naive" ]; then
  echo "running naive BWT implementation on small queries..."
  time ./naive/bwt-naive $BWT_FILE $SMALL_QUERY $PP_DATA $SMALL > $OUTPUT_FILE
  echo "...done running"
  tail -1 $OUTPUT_FILE
fi

if [ "$1" = "run_rmi" ]; then
  echo "running RMI implementation on small queries..."
  t=1
  time ./opt/bwt-match-checkpoint.8bit.32 $BWT_FILE $SMALL_QUERY $PP_DATA 24 $SMALL 0 $t $RMI_OUTPUT_FILE $SMALL > $OUTPUT_FILE
  echo "...done running"
  cat $OUTPUT_FILE
fi

if [ "$1" = "run_lut" ]; then
  echo "running LUT implementation on small queries..."
  t=1
  time ./opt/bwt-match-checkpoint.8bit.32 $BWT_FILE $SMALL_QUERY $PP_DATA 24 $SMALL 0 $t $LUT_OUTPUT_FILE 8 > $OUTPUT_FILE
  echo "...done running"
  cat $OUTPUT_FILE
fi

if [ "$1" = "run_opt" ]; then
   echo "running optimized BWT implementation on small queries..."
   t=4
   time ./opt/bwt-clean-match-checkpoint.8bit.32 $BWT_FILE $SMALL_QUERY $PP_DATA 24 $SMALL 0 $t > $OUTPUT_FILE
   echo "...done running"
   cat $OUTPUT_FILE
fi
