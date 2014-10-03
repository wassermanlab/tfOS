// test_maps.cpp
// g++ -o test  test_maps.cpp
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

/*int main(int argc, char* argv[])
{
	if( argc != 2 ){
        cout << "Missing arguments:" << endl;
        cout << "./test  arg1 " << endl;
        cout << "arg1  /absolute path/INPUT_parameters.txt" << endl;
        exit(1);
    }
	char* infilename = argv[1];  // parameter file to run simulation
	ifstream myinfile(infilename); // set in filestream to the file infilename
*/
int main()
{
	// http://choorucode.com/2011/01/11/c-insert-or-update-a-map/
	// map.insert() will only create, not overwrite a key that exists; so if a key exists insert silently
	//    fails. However insert is more efficient than [] for new entries
	//    I tried the new emplace() but it failed, so I guess we don't have either updated C++ or whatever version is needed
	// if a key exists, using the [] is the best way to update the value assigned to that key
	// if you know a key exists use map.at() to retrieve the value; if it doesnt exist get an out_of_range error 
	// if you check it exists and then do something, use map.find instead. If map.find() is equal to map.end() 
	//      then the key doesnt exist, otherwise find returns an iterator pointing at the element found

	cout << "-----simple map (string + int)"<<endl;
	map< string, int > testmap0;
	testmap0["john"]=  25;
	testmap0["finn"]=  3;
	// insert is 'safe' as will not overwrite a key value pair that already exists
	testmap0.insert( make_pair( "string", 6 ) ); // apparently this way uses type constructors more?
	testmap0.insert(map<string,int>::value_type("string2", 7) ); // apparently this way uses type constructors less?
	cout << "map0: "<< testmap0["john"]<<endl;
	cout << "map0: "<< testmap0["string"]<<"  using insert make_Pair"<<endl;
	cout << "map0: "<< testmap0["string2"]<< "  using insert value_type"<<endl;
	// find is 'safe' as won't create an entry if key doesn't exist. [] will create an entry if key doesn't exist
	if ( testmap0.find("finn") == testmap0.end() ) {
		cout << "finn not found" <<endl;
	}else{
		cout << "finn found" <<endl;
	}
	// OR
	map< string, int >::iterator iter = testmap0.find("finn"); //iter is a pointer
	if(iter==testmap0.end() ){
		cout << "iterator didn't find"<<endl;
	}else{
		cout << "found by iter: "<< iter->first<<" "<<iter->second<<endl;
	}
	cout << "rock ne: "<< testmap0["rock"]<< " end"<<endl; // rock isnt in map, but cout gives '0'
	cout << "use at (should fail): "<< testmap0.at("rock")<< " end"<<endl; // it doesn't fail with out_of_range because apparently rock now exists

	cout <<endl<< "-----simple pair"<<endl;
	pair< string,int > pair1 ("gus" , 10);
	cout <<"pair1: "<< pair1.first<<" " << pair1.second<<endl;
	
	cout <<endl<< "-----map of types string + pair"<<endl;
	map< string, pair< string ,int > > testmap1;
	testmap1["dog"]= make_pair("dog",4);
	pair1=testmap1["dog"];
	cout << "     two ways to get output"<<endl;
	cout << "map1: "<< pair1.first <<" "<< pair1.second<<endl;
	// the below access also works
	cout << "map1: "<< (testmap1["dog"]).first <<" "<< (testmap1["dog"]).second <<endl;

	cout << "-----map of map string, string + int"<<endl;
	map< string, map<string, int> > testmap2;
	testmap2["john"]["natasha"]= 40;
	cout << "map2: "<< testmap2["john"]["natasha"]<<endl;

	cout << "-----map of map string, string + pair"<<endl;
	map< string, map<string, MyPair > > testmap3; // the > > > MUST be spaced apart
	testmap3["mom"]["dad"]=make_pair(70,80);
	testmap3["mom"]["me"]=make_pair(3,7);
	MyPair pair2 (22,23);
	testmap3["a"]["b"]=pair2;
	cout << "map3.1: "<< (testmap3["mom"]["dad"]).first <<" "<< (testmap3["mom"]["dad"]).second <<endl;
	cout << "map3.2: "<< (testmap3["mom"]["me"]).first <<" "<< (testmap3["mom"]["me"]).second <<endl;
	cout << "map3.3: "<< ((testmap3.at("a")).at("b")).first <<" "<< ((testmap3.at("a")).at("b")).second << "   (using  map.at() )" <<endl;

	cout <<endl<< "----- testing insert map of maps"<<endl;
	map<string, MyPair > inner;
	inner.insert(make_pair("stringinner", make_pair(23.4, 11)) );	
	testmap3.insert( map<string, map<string, MyPair > >::value_type( "string3", inner ) ); // apparently this way uses type constructors more?
	cout << "map3.4: "<< ((testmap3.at("string3")).at("stringinner")).first <<" "<< ((testmap3.at("string3")).at("stringinner")).second << "   (using  map.at() )" <<endl;
	
	cout <<endl<< "-----try iterator on map of maps"<<endl;
	MyMapmap::iterator iter1;  //iter is a pointer
	iter1 =  testmap3.find("mom");
	inner.clear();
	if(iter1==testmap3.end() ){
		cout << "mom not found"<<endl;
	}else{
		cout<<"mom found outer map"<< iter1->first<<endl;
		//inner = iter1->second;
		//map<string, MyPair >::iterator iter2 = inner.find("dad"); //this works too
		map<string, MyPair >::iterator iter2  = (iter1->second).find("dad");
		if( iter2 !=  (iter1->second).end() ){
			// the inner map is string + pair, so iter2 needs 2nd element of map but both first and 2nd element of pair
			cout << "'dad' found inner map: "<< (iter2->second).first <<" "<< (iter2->second).second <<endl;
		}else{
			cout << "iter2 did not find 'dad'" <<endl;
		} 
	}

	cout <<endl<< "-----typedef of map of map string, string + pair"<<endl;
	MyMapmap testmap4; // typedef above
	testmap4["paul"]["jutta"]=make_pair(40,44);
	cout << "map4: "<< (testmap4["paul"]["jutta"]).first <<" "<< (testmap4["paul"]["jutta"]).second <<endl;
	
	cout << endl << "-----test if key exists in testmap3 "<<endl;
	if( testmap3.count("mom") >0){
		cout << "'mom' exists "<< testmap3.count("mom") <<endl;
	}
	cout << "'dad' not in outer map "<< testmap3.count("dad") <<endl;
	cout << "'dad' in inner map "<< (testmap3["mom"]).count("dad") <<endl;


	return 0; 
}
