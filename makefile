all:
	g++ main.cpp -o main && ./main src/test_genome.txt src/sample_hw_dataset.txt
clean:	
	rm main

