#!/usr/bin/env python
# coding: utf-8
from sys import path
import urllib.request
from bs4 import BeautifulSoup
import time 
import pandas as pd
import argparse
import os


RUN = "SRR10211551"
FORMATS = [".bam",".fastq",".fq"]
AWS = ["s3"]

def get_url(run,formats=FORMATS):

    page = urllib.request.urlopen('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run='+run)
    p = page.read()

    soup = BeautifulSoup(p)

    aws_address = []
    for item in soup.find_all('tbody'):
        
        # The tbody have AWS info 
        if((str(item).find("AWS")>=0)):
            
            # Url is stored in <a href=>
            a_url = (item.find_all('a', href=True))
            
            # Found link in the <a href>
            if(len(a_url)>0):
                for a in a_url: 
                    if (any((str(a).find(f)>=0) for f in formats) & \
                        (str(a).find("s3")>=0)): 
                            url = a['href']
                            print(url)
                            aws_address.append(url)

            # Url is stored in <td>
            else:
                item_td = item.find_all('td')
                for i in item_td:
                    if (any((str(i).find(f)>=0) for f in formats) & \
                        (str(i).find("s3")>=0)): 
                            url = (str(i).strip('<td>').strip("</"))
                            print(url)
                            aws_address.append(url)

    s3_address=[]
    https_address=[]

    for url in aws_address:
        if((url).find("https://")>=0): 
            https_address.append(url)
            s3_url = url.replace("https://", "s3://")
            s3_url = s3_url.replace(".s3.amazonaws.com", "")
            
            print(s3_url)
            s3_address.append(s3_url)
        elif((url).find("s3://")>=0):
            s3_address.append(url)
            https_url = url.replace("s3://","https://")
            https_url = https_url.replace("/SRR", ".s3.amazonaws.com/SRR")
            https_address.append(s3_url)


    df_result = pd.DataFrame({"s3":s3_address,"http":https_address})
    df_result['run'] = run

    return df_result

def run_main(args):

    runs = args.runs

    df_result = []
    for r in runs:
        time.sleep(args.sleep)
        df = get_url(r,formats=FORMATS)

        if(len(df_result)<1):
            df_result= df
        else:
            df_result=df_result.append(df)

    df_result.to_csv(os.path.join(args.outpath,"query_result.csv"),index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert SRR to a set of down load address.')
    parser.add_argument('--runs', nargs='+',default=["SRR10211551","SRR10211552"], type=str)
    parser.add_argument('--outpath', '-o',default="./",  type=str,  help='An output folder')
    parser.add_argument('--sleep', '-s',  type=int,default=2,  help='Sleeping time for each query')

    args, unknown = parser.parse_known_args()
    run_main(args)