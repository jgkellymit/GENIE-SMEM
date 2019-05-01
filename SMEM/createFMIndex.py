from os import path
import json


def createFMIndex(referenceSequenceFile):
    with open(path.join("data", referenceSequenceFile), "r") as refFile:
        refFile.readline()
        og_sequence = ""
        for line in refFile:
            og_sequence += line.strip()

    bwt_array, first_index_array, suffix_array = createBWTMatrix(og_sequence)

    o_matrix = getOccuranceMatrix(bwt_array)
    c_dic = getCountDic(first_index_array)
    with open(path.join("data", referenceSequenceFile.split(".")[0] + "-FM.json"), "w") as fm_file:
        fm_file.write(json.dumps({"bwt_array":bwt_array, "suffix_array": suffix_array, "occurance_matrix": o_matrix,
                                  "count_dic": c_dic}, indent=4, sort_keys=True))


def createBWTMatrix(ref_seq):
    ref_seq += "$"
    ref_size = len(ref_seq)
    bwt_matrix = []

    # Rotate
    for index in range(ref_size):
        bwt_matrix.append(ref_seq[index:] + ref_seq[:index])
    bwt_matrix.sort()

    bwt_array = []
    first_index_array = []
    suffix_array = []
    for rotation in bwt_matrix:
        bwt_array.append(rotation[-1])
        first_index_array.append(rotation[0])
        suffix_array.append(ref_size - rotation.index("$"))

    return bwt_array, first_index_array, suffix_array


def getOccuranceMatrix(BWT_array):
    o_dic = {}
    index = 0
    for item in BWT_array:
        if item not in o_dic:
            o_dic[item] = [0]
        item_list = o_dic[item]
        while len(item_list) <= index:
            item_list.append(item_list[-1])
        item_list[-1] += 1

        index += 1

    # Extend everything to full matrix length
    for key in o_dic:
        item_list = o_dic[key]
        while len(item_list) < index:
            item_list.append(item_list[-1])

    return o_dic


def getCountDic(first_index_array):
    count_dic = {}
    index = 0
    for char in first_index_array:
        if char not in count_dic:
            count_dic[char] = index
        index += 1
    count_dic[""] = index
    return count_dic


createFMIndex("mississippi.fa")

if __name__ == "__main__":
    import sys

    filename = "mississippi.fa"

    if len(sys.argv) > 1:
        filename = sys.argv[1]

    createFMIndex(filename)