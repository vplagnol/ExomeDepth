#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_overlap (const SEXP chromosome, const SEXP start, const SEXP end, const SEXP index_chromosome, const SEXP index_start, const SEXP index_end, const SEXP index_name, const SEXP fraction_overlap_required); 
}




class indel_call {
public:
  int start, stop, length, count, id;
  string chrom, type;
  string line;  //the full line that describes this feature

  indel_call (const string chrom_a, const int start_a, const int stop_a, const string & line_a, const int id_a): chrom (chrom_a), start (start_a), stop (stop_a),line (line_a), id(id_a)  {};
  
  bool operator< (const indel_call &compare) const;
};



bool indel_call::operator<(const indel_call& other) const {
  if (chrom < other.chrom) {return true;}
  if (chrom > other.chrom) {return false;}
  if (start < other.start) {return true;}
  if (start > other.start) {return false;}
  if (type < other.type) {return true;}  //now same position, same chromosome, then if type is smaller then indeed it will be smaller
  return (id < other.id);
  //return false;  
}


SEXP C_overlap (const SEXP chromosome, const SEXP start, const SEXP end, const SEXP index_chromosome, const SEXP index_start, const SEXP index_end, const SEXP index_name, const SEXP fraction_overlap_required) {

  const int nqueries = 100;
  
  SEXP final;
  PROTECT(final = NEW_CHARACTER(nqueries)); 
  UNPROTECT(1);

  return (final);
}







/*
int main (int narg, char ** argc) {
  clock_t c_start = clock();

  cout<<"Running overlap, use ./overlap -h for more info\n";
 
  parameters * my_params = new parameters(narg, argc);


  vector<set< indel_call > > list_indels (255);


  cout<<"Parsing file "<<my_params->input_file<<endl;  
  ifstream inp (my_params->input_file.c_str());
  if (!inp) {cerr<<"Cannot open "<<my_params->input_file.c_str()<<endl;exit(1);}
  
  //------------------------------------------------------------ reads the index of features
  istringstream instream;
  string line;
  string chrom, gene, id;
  int start, end;
  int nindex = my_params->geneIndex.size();

  vector< map<string, vector<int> > > map_start (nindex), map_end(nindex);
  vector< map<string, vector<string> > >  map_genes (nindex);
  
  
  for (int index = 0; index != my_params->geneIndex.size(); index++) {
    
    cout<<"Index "<<index<<" :"<<my_params->geneIndex[index]<<endl;
    ifstream ind(my_params->geneIndex[index].c_str());
    if (!ind) {cerr<<"Cannot open gene index "<<my_params->geneIndex[ index ]<<endl;exit(1);}
    
    while (1) {
      getline(ind, line);
      instream.clear(); 
      instream.str(line);
      if (ind.eof()) break;
      instream>>chrom>>start>>end>>id;
      map_start[index][ chrom ].push_back(start);
      map_end[index][ chrom ].push_back(end);
      map_genes[index][ chrom ].push_back(id);
    }

    ind.close();

     ///----------- now sorting, need to do each for each chromosome

    for (map<string, vector<int> >::iterator it = map_start[index].begin(); it != map_start[index].end(); ++it) {
     string chrom = (*it).first;

     int nelems =  map_end[index][ chrom ].size();
     vector<int> indexes (nelems, 0);
     for (int i = 0; i != nelems; i++) indexes[i] = i;

     sort( indexes.begin(), indexes.end(), index_cmp<vector<int>&>( map_start[index][ chrom  ]));
     vector<int> temp1(nelems, 0);
     vector<int> temp2(nelems, 0);
     vector<string> temp3(nelems);

     for (int i = 0; i != nelems; i++)   {
        temp1[ i ] =  map_start[index][ chrom ][ indexes[i] ];
        temp2[ i ] =  map_end[index][ chrom ][ indexes[i] ];
        temp3[ i ] =  map_genes[index][ chrom ][ indexes[i] ];
      }

     (map_start[index][ chrom ]) = temp1;
     (map_end[index][ chrom ]) = temp2;
     (map_genes[index][ chrom ]) = temp3;
    }

  }
  cout<<"Done reading the indexes\n";


  //////////------------------------------------------------------- reads the features we care about

  string Chrom, Start, Stop, headers;
  int ncalls = 0;
  while (1) {
    
    getline(inp, line);
    instream.clear(); 
    instream.str(line);
    
    if (inp.eof()) break;
    if (!my_params->onebp) instream>>Chrom>>Start>>Stop; 
    if (my_params->onebp) instream>>Chrom>>Start;
    
    if ( (Chrom == "#CHROM") || (Chrom == "Chr") || (Chrom == "chromosome") || (Chrom == "seqname") || (Chrom == "ReadName")) {
      headers = line;
    } else {
      if (Chrom[0] == '#') continue;

      if (Chrom.substr(0, 3) == "chr") Chrom = Chrom.substr(3);

      ncalls++;
      int my_tag = 0;
      int my_start = atoi(Start.c_str());
      int my_stop;
      if (!my_params->onebp) my_stop = atoi(Stop.c_str());
      if ( my_params->onebp) my_stop = my_start + 1;
      
      indel_call * my_call = new indel_call (Chrom, my_start, my_stop, line, ncalls);
      list_indels[ my_tag ].insert(*my_call);
	//cout<<list_indels[0].size()<<"  "<<my_tag<<endl;
    }
  }
  
  cout<<"Number of calls: "<<ncalls<<", "<<list_indels[0].size()<<" of them being unique"<<endl;
  inp.close();

  //----------------------------------------------------------------- Now output everything
  cout<<"Placing output in "<<my_params->output_file<<endl;
  ofstream out (my_params->output_file.c_str());

  ///---- headers
  out<<headers;
  for (int index = 0; index != nindex; index++) {
    string myname = my_params->geneIndex[ index ];
    unsigned int pos =  myname.rfind ('/');
    out<<"\t"<<myname.substr( pos + 1, myname.size() - pos - 1);
    //out<<"\t"<<my_params->geneIndex[ index ];
  }
  out<<"\n";
  //----


  for (int tag = 0; tag != 1; tag++) {
    for (set<indel_call>::iterator i = list_indels[tag].begin(); i != list_indels[tag].end(); ++i) {  //for each deletion
      
      //out<<(*i).chrom<<"\t"<<(*i).start<<"\t"<<(*i).stop<<"\t";
      out<<(*i).line;

      //-------------------------------------
      for (int index = 0; index != nindex; index++) {  //for each index
	
	string closeGene="NA";		  
	if (map_start[index].count( (*i).chrom ) > 0) {
	  
	  int pos = upper_bound ( map_start[index][ (*i).chrom ].begin(), map_start[index][ (*i).chrom ].end(), (*i).start) -  map_start[index][ (*i).chrom ].begin();
	  int wide = 100;
	  pos = pos - wide;
	  int maxRange = pos + 2*wide;
	
	  if (pos < 0) pos = 0;
	  if (maxRange > map_start[index][ (*i).chrom ].size()) maxRange = map_start[index][ (*i).chrom ].size();      
	  bool flag = false;
	  
	  int mypos = pos;
	  for (int mypos = pos; mypos != maxRange; mypos++) {
	    
	    int overlap = max(0, min(map_end[index][ (*i).chrom ][ mypos ], (*i).stop) - max(  map_start[index][ (*i).chrom ][ mypos ], (*i).start ));
	    int threshold = 1;

	    if (( !my_params->onebp) && ( my_params->min_overlap[ index ] > 0)) {
	      //threshold = (int) max( my_params->min_overlap * ((double) ((*i).stop - (*i).start)),  my_params->min_overlap* ( (double) (map_end[index][ (*i).chrom ][ mypos ] - map_start[index][ (*i).chrom ][ mypos ])));
	      threshold = (int) ( my_params->min_overlap[ index ] * ((double) ((*i).stop - (*i).start)));  //the ref feature must be at least a percentage of the called feature
	      //threshold = (int) (my_params->min_overlap * ( (double) (map_end[index][ (*i).chrom ][ mypos ] - map_start[index][ (*i).chrom ][ mypos ])));
	    }
	    
	    //cout<<"--------------- "<<overlap <<"   "<< threshold<<" "<<map_end[index][ (*i).chrom ][ mypos ]<<"\t"<<endl;

	    if ( overlap >= threshold ) {
	      if (!flag) closeGene = map_genes[index][ (*i).chrom ][ mypos ]; else closeGene += "," + map_genes[index][ (*i).chrom ][ mypos ];	      
	      flag = true;	      
	      maxRange++;
	      if (maxRange > map_start[index][ (*i).chrom ].size()) maxRange = map_start[index][ (*i).chrom ].size();      
	      //cout<<flag<<"\t"<<mypos<<"\t"<<maxRange<<endl;
	    }
	    if (flag == false) closeGene = "none";
	    //if ( map_start[index][ (*i).chrom ][ mypos ] > (*i).stop ) break;
	  }
	} else {closeGene = "none";}
	out<<"\t"<<closeGene;
      }
      out<<endl;
      //------------------------------------- 
    }

  }
  out.close();



  clock_t finish = clock();
  cout<<"Time needed, in seconds: "<<(double(finish)-double(c_start))/CLOCKS_PER_SEC<<endl;


}
*/
