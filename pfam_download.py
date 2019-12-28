import requests
from bs4 import BeautifulSoup
import os
from os import path
import re
import socket
from time import gmtime, strftime
import sys

def is_connected():
    try:
        socket.create_connection(("www.google.com", 80))
        return True
    except OSError:
        pass
    return False

def accession_list():

	browse_list = 'abcdefghijklmnopqrstuvwxyz'
	browse_list = list(browse_list) + ['numbers', 'new', 'top%20twenty']
	line_list = []
	accession_with_size = {}
	total = 0

	for key in browse_list:
		print('Fetching entries that start with: ' + key)
		browse_url = 'https://pfam.xfam.org/family/browse?browse=' + key
		response = requests.get(browse_url)

		if response.status_code == 200:
			soup = BeautifulSoup(response.text, 'html.parser')
			find = soup.find_all('tr')
			for i in find:
				line = i.text
				line = line.split('\n')
				for l in line:
					if l != '':
						line_list.append(l)

		for i in range(len(line_list)):
			x = line_list[i]
			if re.search('^PF', x):
				try:
					size = int(line_list[i + 3])
					if size >= 400:
						if x in accession_with_size:
							continue
						else:
							total += 1
							accession_with_size[x] = size
				except:
					continue

	with open('resources/accession_list.txt', 'w') as f:
		for i in accession_with_size:
			f.write(str(i) +"\n")

	print('pickle dumped......')

def download_entries(accession_with_size):
	for key in accession_with_size:
		if is_connected():
			filename = 'families/' + key + '.fasta'
			if path.exists(filename) is False:
				download_url = 'https://pfam.xfam.org/family/' + key + '/alignment/full/format?format=fasta&alnType=full&order=t&case=l&gaps=dashes&download=0'
				try:
					print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
					print('Downloading alignment: ' + key)
					response = requests.get(download_url)
					soup = BeautifulSoup(response.text, 'html.parser')
					text = soup.get_text()
					filename = 'families/' + key + '.fasta'
					with open(filename, 'w') as f:
						for line in text:
							f.write(line)
				except:
					continue
		else:
			print('Internet disconnected!')
			sys.exit()

def main():
	#accession_list()
	#sys.exit()
	accession_with_size = []
	with open('resources/accession_list.txt', 'r') as f:
		for line in f:
			accession_with_size.append(line.strip('\n'))
	#print(accession_with_size, len(accession_with_size))
	download_entries(accession_with_size)

if __name__ == '__main__':
    main()
