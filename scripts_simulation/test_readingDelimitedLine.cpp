// test_readingDelimitedLine.cpp
// g++ -o test  test_readingDelimitedLine.cpp
// ./test  filename

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	if( argc != 2 ){
        cout << "Missing arguments:" << endl;
        cout << "./test  arg1 " << endl;
        cout << "arg1  /absolute path/INPUT_parameters.txt" << endl;
        exit(1);
    }
    char* infilename = argv[1];  // parameter file to run simulation
	ifstream myinfile(infilename); // set in filestream to the file infilename
	char cline[101];
	//string line;  // string of each line            
	int flag=0;
	//while( getline(myinfile, line))
	while (!myinfile.eof())
	{
		myinfile.getline(cline, 101);
		string line(cline);
		/*string line;
		getline(myinfile, line); */ // put each fileline into line 
		istringstream ss(line);  //line becomes an input stream; Then each direction operator >> directs the next element in the stream
								// can then get each word in line; as character based streams are deemed to be separated by whitespace \n,\t,' '
		string first, second;
		int var1, var2;
		int size=0;
		// counting the number of words, but it doesn't reset the ss iterator? back to start of ss
		/*string word;
		while( ss >>word ){
			cout <<"word "<<word <<endl;
			size++;
		}
		cout << "size "<< size<<endl;
		*/

		ss >> first >> second >> var1 >> var2; // here I know there are 4 words per line
		if( first.compare( "coop_interactions") == 0 || line.empty()){ 
			flag=1;
			continue;
		}
		if( first.compare("end_interactions")== 0){ 
			break;
		}
		if(flag==1){
			cout <<"first "<< first ;
			cout <<" second "<< second ;
			cout << " third "<< var1;
			cout << " fourth "<< var2<<endl;
		}
		

		cout << endl;
		/*// if I don't know how many lines and they are all the same type, can loop over them
		  //   I also explored converting a string to a double 
		//while( !ss.eof()){ // this allows blank lines in file to be in ss
		int count=0;
		string str;
		while( ss >> str ){ // this doesnt assign blank lines
			//ss >> first;  // used with while(!ss.eof())
			// now could assign each element to a vector or something
			count++;
			if(count<4) {
				cout << "each ss "<< str <<endl;
			}else{
				count=0;
				double num;
				ss >> str;
				num = atof(str.c_str());
				cout << "string "<< str << " num "<<num<< endl;
			}
		}*/

	}
	myinfile.close();
 return 0; 
}
