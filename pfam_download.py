import requests
from bs4 import BeautifulSoup
import os
from os import path
import re
import socket
from time import gmtime, strftime
import sys
import subprocess

#function pings google.com to check if the computer is connected to the internet
def is_connected():
    try:
        socket.create_connection(("www.google.com", 80))
        return True
    except OSError:
        pass
    return False
#function for creating a list of accession number from pfam's entry table
def accession_list():
        #pfam entry table is organized alphabetically
        browse_list = 'abcdefghijklmnopqrstuvwxyz'
        #pfam entry table also has these three fields apart from alphabets
        browse_list = list(browse_list) + ['numbers', 'new', 'top%20twenty']
        line_list = []
        accession_with_size = {}
        total = 0
        

        for key in browse_list:
                print('Fetching entries that start with: ' + key)
                browse_url = 'https://pfam.xfam.org/family/browse?browse=' + key
                #a response object for the above url
                response = requests.get(browse_url)
                #if there is a good repsonse
                if response.status_code == 200 or response.status_code == '200':
                        try:
                                #BeautifulSoup parses the html part of the text from the response
                                soup = BeautifulSoup(response.text, 'html.parser')
                                #for finding all cells with a tr tag in the entry table 
                                find = soup.find_all('tr')
                                for i in find:
                                        line = i.text
                                        line = line.split('\n')
                                        for l in line:
                                                #if row is not empty then append it to line_list
                                                if l != '':
                                                        line_list.append(l)
                        except:
                                print('Response Error')
                                sys.exit()
                else:
                        continue
                
                #for each row in line_list
                for i in range(len(line_list)):
                        x = line_list[i]
                        #if row has an accession number then proceed 
                        if re.search('^PF', x):
                                try:    
                                        #size of the family is after every third element after accession
                                        size = int(line_list[i + 3])
                                        #get families only with more than 500 sequences
                                        if size >= 500:
                                                #if family is already in list then continue
                                                if x in accession_with_size:
                                                        continue
                                                else:
                                                        #if family is not in list then add it to the dictionary
                                                        total += 1
                                                        accession_with_size[x] = size
                                except:
                                        continue

        #write all accession numbers in a file
        with open('temp_files/accession_list.txt', 'w') as f:
                for i in accession_with_size:
                        f.write(str(i) +"\n")

        print('pickle dumped......')

#downloading all entries in the accession list file
def download_entries(accession_with_size):
        for key in accession_with_size:
                if is_connected():
                        filename = 'pfam_entries/' + key + '.fasta'
                        if path.exists(filename) is False:
                                #this line is from when I restarted the download to another folder (from families folder to pfam_entries folder)
                                old_filename = 'families/' + key + '.fasta'
                                if path.exists:
                                        #delete old file before downloading new file
                                        os.system('rm -rf ' + old_filename)
                                download_url = 'https://pfam.xfam.org/family/' + key + '/alignment/full/format?format=fasta&alnType=full&order=t&case=u&gaps=dashes&download=0'
                                try:    
                                        #print current date and time
                                        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                                        print('Downloading alignment: ' + key)
                                        response = requests.get(download_url)
                                        if response.status_code == 200 or response.status_code == '200':

                                                soup = BeautifulSoup(response.text, 'html.parser')
                                                text = soup.get_text()
                                                filename = 'pfam_entries/' + key + '.fasta'
                                                with open(filename, 'w') as f:
                                                        for line in text:
                                                                if line.startswith('500 Internal Server Error'):
                                                                        return 
                                                                f.write(line)
                                        else:
                                                continue
                                except:
                                        continue
                else:
                        return

def main():

        #print('If you think that you have enough families downloaded then just press Ctrl + Z and stop the script :)')
        accession_list()
        accession_with_size = []
        with open('temp_files/accession_list.txt', 'r') as f:
                for line in f:
                    #appending all accessions numbers in a list
                    accession_with_size.append(line.strip('\n'))

        '''
        dict_size = len(accession_with_size)
        while True:
                download_entries(accession_with_size)
                os.chdir('media/Data/consensus/pfam_entries')
                size_of_directory = int(subprocess.check_output('ls -l | wc -l', shell=True))
                if dict_size == size_of_directory:
                        break
        print('Download Complete!')
        sys.exit()
        '''

if __name__ == '__main__':
    main()
