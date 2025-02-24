import json
import sys
import os

def summarize(folder):
    summary_dict = {}

    for filename in os.listdir(folder):
        if filename != "summary.json": 

            file_path = os.path.join(folder, filename)
            
            with open(file_path, 'r') as json_file:
                subset_dict = json.load(json_file)

            summary_dict[os.path.splitext(filename)[0]] = subset_dict

    summary_json = json.dumps(summary_dict, indent=4)
    file_name = "summary.json"
    file_path = os.path.join(folder, file_name)

    with open(file_path, 'w') as json_file:
        json_file.write(summary_json)       
    print("finished summarizing: ", folder)

if __name__ == '__main__':
    folder = sys.argv[1]
    summarize(folder)