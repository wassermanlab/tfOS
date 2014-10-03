
// test_char.cpp
// g++ -o test  test_char.cpp
// ./test 

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
//#include <unordered_map>
#include <map>
#include <utility>  //for pair and make_pair
#include <vector>

using namespace std;
typedef std::pair< float,int> MyPair;
typedef std::map<string, map<string, MyPair > > MyMapmap; // the > > > must be spaced apart

int main(int argc, char* argv[])
{
	if( argc != 2 ){
        cout << "Missing arguments:" << endl;
        cout << "./test  arg1 " << endl;
        cout << "arg1  inputfile" << endl;
        exit(1);
    }
	char* infilename = argv[1];  // parameter file to run simulation
	const int BUFF=100;
	char buf[BUFF];
	char* bufpt;
	ifstream myinfile(infilename); // set in filestream to the file infilename
 	if( myinfile.is_open()) {
		cout << "file is open: " <<infilename<<endl;
		while ( !myinfile.eof()) {
			myinfile.getline(buf, BUFF);
			cout << "line retrieved: "<<" "<<buf <<endl;
			bufpt = buf;
			if( strlen(bufpt)>5){
				cout << "length "<<strlen(bufpt) << endl;
			}else{
				cout << "too short ("<< strlen(bufpt)<<")"<<endl;
			}
		}
	}
	string len;
	if(len==NULL){
		cout << "len is null"<<endl;
	}
	
	/*
	cout << "test how strlen() handles" << endl;
	string teststr;
	teststr = "12characters";
	cout << "len "<< teststr.size()<<endl;
	char testchar[] ="12characters";
	char* tmp = testchar;
	//char* tmp = "12characters";
	cout << "char len "<<strlen(tmp) <<endl<<endl; // strlen needs char*
	*/
	/*
	cout << "test string conversion to double<<endl;
	string strnum = "19.9";
	double num;
	num = atof(strnum.c_str());
	cout << "str "<< strnum << " num "<< num<<endl;
	*/

	return 0;
}
