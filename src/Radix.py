radix = 0;
num = 0;
str = 'GGNC'
n = len(str)
for i in range(len(str)):
	if str[i] == 'A':
		num = 0
	elif str[i] == 'C':
		num = 1
	elif str[i] == 'G':
		num = 2
	elif str[i] == 'T':
		num = 3
	radix = radix + num * pow(4, n-1)
	n = n -1;
print(radix)
