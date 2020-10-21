import json
import requests, sys
import click

def print_coords_from_json(obj):
    chromosome = "chr" + str(obj['seq_region_name'])
    if 'Transcript' in obj.keys():
        for trs in obj['Transcript']:
            display_name = trs['display_name']
            file_suffix = display_name+"_"+trs['id']+'.bed'
            print("writing ",display_name,trs['id'])
            if trs['is_canonical'] > 0: # canonical transcript
                # if canonical, insert C_ prefix
                fh = open("C_"+file_suffix, 'w')
            else:
                fh = open(file_suffix, 'w')
            for exon in trs['Exon']:
                fh.write(chromosome +
                         "\t" + str(exon['start']) +
                         "\t" + str(exon['end']) +
                         "\t" + display_name +
                         "\n")
    else:
        # it is a single transcript only
        for exon in obj['Exon']:
            intree.add(Interval(exon['start'], exon['end']))

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--id', '-i', type=str, help='ENSEMBL ID to print out', required=True)
@click.option('--rest/--no-rest', '-r', type=bool, help='Read from ENSEMBL REST or JSON file',
              required=False, default=True)
def printJSON(id, rest):
    obj = None  # the JSON object we are playing with
    if rest:
        server = "https://rest.ensembl.org"
        ext = "/lookup/id/" + id + "?expand=1"
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        obj = r.json()
        json_file = open(id + '.json', 'w')
        json_file.write(r.content.decode('utf-8'))
        json_file.close()
    else:   # we are reading from file
        with open(id + '.json', 'r') as myfile:
            data = myfile.read()
        obj = json.loads(data)
    print_coords_from_json(obj)

if __name__ == "__main__":
    printJSON()