
import requests
import sys

ORG = sys.argv[1]

url = 'https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa'
params = {
    'boolshad.sendToGalaxy': '0',
    'boolshad.sendToGenomeSpace': '0',
    'boolshad.sendToGreat': '0',
    'clade': 'mammal',
    'db': '0',
    'hgsid': '611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa',
    'hgta_compressType': 'none',
    'hgta_doTopSubmit': 'get output',
    'hgta_group': 'allTracks',
    'hgta_outFileName': 'rmsk.gtf',
    'hgta_outputType': 'gff',
    'hgta_regionType': 'genome',
    'hgta_table': 'rmsk',
    'hgta_track': 'rmsk',
    'jsh_pageVertPos': '0',
    'org': ORG,
    'position': ''}
       
with requests.Session() as session:
    response = session.post(url, data=params) 
        
sys.stdout.write(response.content.decode('utf-8'))
