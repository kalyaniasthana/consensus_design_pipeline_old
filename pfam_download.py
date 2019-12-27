import requests
from bs4 import BeautifulSoup 


'''
fasta_url = 'https://pfam.xfam.org/family/PF10417/alignment/full/format?format=fasta&alnType=full&order=t&case=l&gaps=dashes&download=0'

response = requests.get(fasta_url)
if response.status_code == 200:
	print('you are good to go')
	soup = BeautifulSoup(response.text,'html.parser')
	text = soup.get_text()
	with open('test_file.fasta', 'w') as f:
		for line in text:
			f.write(line)

'''

browse_list = 'abcdefghijklmnopqrstuvwxyz'
browse_list = list(browse_list) + ['numbers', 'new', 'top%20twenty']
line_list = []
accession_with_size = {}
total = 0

for key in browse_list:
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
		if x.startswith('PF'):
			try:
				size = int(line_list[i + 3])
				if size >= 500:
					if x in accession_with_size:
						print(key, x, 'ALREADY THERE')
					else:
						print(key, x, size)
						total += 1
						accession_with_size[x] = size
			except:
				continue

print(len(accession_with_size), 'DICT********************')
print(total, 'TOTAL########################')