#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import sys


def get_sample_info_from_json(record_name, json_filename, runID=None):

    sampleinfo = {"json_file": None, "episode": None, "accession" : None, "run": None, "lib_repeat" : None}

    path, file = os.path.split(json_filename)

    sampleinfo["json_file"] = file
    sampleinfo["run"] = runID
    
    dataset_id = re.sub('_S\\d*_L\\d*$', '', record_name)

    dataset_id = re.sub('WCMID-', '', dataset_id)

    for identifier in dataset_id.split("-"):
        if identifier.startswith("A"):
            sampleinfo["accession"] = identifier
        elif identifier.startswith("C") or "P" in identifier or "NCTC" in identifier:
            sampleinfo["accession"] = "POSCONTROL"
        elif identifier[0].isdigit:
            epi_repeat = identifier.split("R")
            if len(epi_repeat) > 1: 
                sampleinfo["episode"] = epi_repeat[0]
                sampleinfo["lib_repeat"] = epi_repeat[1]
            else:
                sampleinfo["episode"] = epi_repeat[0]

    return  {f'1sample.{k}': v for k, v in sampleinfo.items()}

def read_json(json_file):
    with open(json_file) as f:
        json_data = json.load(f)

    return json_data

def parse_phylo_data_from_json(json_data):
    phylo_dict = {}

    phylo_data = json_data.get("phylogenetics")

    for prop in list(phylo_data):
        s_phylolist = sorted(phylo_data.get(prop), key=lambda x: (phylo_data.get(prop)[x]['percent_coverage']), reverse=True)
        for val, list_prop in enumerate(s_phylolist, 1):
            prop_val = prop + str(val)
            prop_pct = prop_val + "_pct" 
            prop_med = prop_val + "_median" 

            phylo_dict[prop_val] = list_prop
            phylo_dict[prop_pct] = phylo_data.get(prop).get(list_prop).get("percent_coverage")
            phylo_dict[prop_med] = phylo_data.get(prop).get(list_prop).get("median_depth")

    return {f'2phylo.{k}': v for k, v in phylo_dict.items()}

def parse_amr_data_from_json(json_data):
    amr_dict = {}

    amr_data = json_data.get("susceptibility")

    for drug in list(amr_data):
        if "S" in amr_data.get(drug).get("predict"):
            amr_dict[drug] = amr_data.get(drug).get("predict")
        else:
            variants = list(amr_data.get(drug).get("called_by"))

            variant_string_list = []

            for variant in variants:
                depth_stats = {}
                mykrobe_filter = amr_data.get(drug).get("called_by").get(variant).get("info").get("filter")
                coverage = amr_data.get(drug).get("called_by").get(variant).get("info").get("coverage")
                depth_stats["ref_depth"] = coverage.get("reference").get("median_depth")
                depth_stats["alt_depth"] = coverage.get("alternate").get("median_depth")
                depth_stats["alt_cov_pct"] = coverage.get("alternate").get("percent_coverage")

                variant_string = variant+ ":" + ":".join(str(x) for x in list(depth_stats.values()))

                variant_string_list.append(variant_string)

            amr_dict[drug] = "R"
            amr_dict[drug + "_determinants(variant:ref_depth:alt_depth:alt_cov_pct)"] = "|".join(variant_string_list)

    return {f'3res.{k}': v for k, v in amr_dict.items()}


def write_data_to_csv(sample_data, output_file):
    header = []

    for record in sample_data:
        for item in list(record):
            if item not in header:
                header.append(item)

    header.sort()

    with open(output_file, 'w') as csv_out:
        dict_writer = csv.DictWriter(csv_out, header)
        dict_writer.writeheader()
        dict_writer.writerows(sample_data)


def parse_args():
    parser = argparse.ArgumentParser(description="Make nice CSV from mykrobe output")
    parser.add_argument("-o", dest="output_file", type=str, help="Path to output CSV file", required=True)
    parser.add_argument("-j", dest="json_files", nargs='+', type=str, help="Path to JSON files", required=True)
    parser.add_argument("-r", dest="sequencing_run_id", type=str, help="Sequencing Run ID", required=True)
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_args()
    
    csv_data = []

    for file in args.json_files:
        json_data = read_json(file)

        for record in json_data:
            all_sample_data = {}
            all_sample_data.update(get_sample_info_from_json(record, file, args.sequencing_run_id))
      
            all_sample_data.update(parse_phylo_data_from_json(json_data.get(record)))
            all_sample_data.update(parse_amr_data_from_json(json_data.get(record)))


        csv_data.append(all_sample_data)
    
    write_data_to_csv(csv_data, args.output_file)
