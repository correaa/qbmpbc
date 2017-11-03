#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#include "../qbox/vector.hpp"
#include "../qbox/cell.hpp"
#include<map>
#include<string>
#include<sstream>
#include<iostream>
namespace qbox{
	using std::string;
	using std::pair;
	using boost::filesystem::path;
	using boost::filesystem::ifstream;
	using boost::filesystem::ofstream;
	using std::clog;
	using std::endl;
	using std::map;
	struct atoms_list : map<string, pair<string, vector> >{
		friend 
		std::ostream& operator<<(std::ostream& os, atoms_list const& self){
			for(atoms_list::const_iterator it=self.begin(); it!=self.end(); ++it)
				os<<"atom "<<it->first<<" "<<it->second.first<<" "<<it->second.second<<endl;
			return os;
		}
	};
	struct species_list : map<string, path>{
		friend 
		std::ostream& operator<<(std::ostream& os, species_list const& self){
			for(species_list::const_iterator it=self.begin(); it!=self.end(); ++it)
				os<<"species "<<it->first<<" "<<it->second<<endl;
			return os;
		}
	};
	class sys{
		species_list species_;
		public:
		atoms_list atoms;
		cell cell_;
		sys(path const& p){
			ifstream ofs(p);
			for(string line; getline(ofs, line);){
				std::istringstream line_ss(line);
				string cmd; line_ss>>cmd;
				if(cmd=="set"){
					line_ss>>cmd;
					if(cmd=="cell"){
						vector a1, a2, a3;
						//double a11, a12, a13, a21, a22, a23, a31, a32, a33;
						line_ss>>a1>>a2>>a3;
						cell_=qbox::cell(a1, a2, a3);
					}
				}else if(cmd=="species"){
					string name, uri; line_ss>>name>>uri;
					species_[name]=path(uri);
				}else if(cmd=="atom"){
					string name, species; vector position; line_ss>>name>>species>>position;
					atoms[name]=std::make_pair(species, position);
				}else{
					clog<<"warning: command "<<cmd<<" not handled in "<<p<<endl;
				}
			}
		}
		sys& operator*=(double factor){
			cell_=qbox::cell(cell_.a(0)*factor,cell_.a(1)*factor, cell_.a(2)*factor);
			for(atoms_list::iterator it=atoms.begin(); it!=atoms.end(); ++it)
				it->second.second*=factor;
			return *this;
		}
		friend
		path const& operator<(path const& p, sys const& self){
			ofstream ofs(p);
			ofs
				<<"set cell ";
					ofs<<self.cell_.a(0)<<endl
				<<self.species_
				<<self.atoms;
			return p;
		}
	};
}//ns qbox
