#ifdef compile_instructions 
ln -sf $0 .$0.cpp && \
c++ -std=c++0x -D_SIESTA_HPP_TEST -I$HOME/prj `#-Wfatal-errors` -Wall -Wno-unused-variable `pkg-config --libs gsl` -lboost_filesystem .$0.cpp -o .$0.x -Wall && \
./.$0.x $@; exit
#endif
#ifndef SIESTA_HPP_
#define SIESTA_HPP_
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#include<boost/spirit/include/qi.hpp>
#include<boost/spirit/include/phoenix_core.hpp>
#include<boost/spirit/include/phoenix_operator.hpp>
#include<boost/iomanip.hpp> //for ignore
#include "boost/units/systems/atomic.hpp"
#include<iomanip>
#include<map>
#include<boost/algorithm/string.hpp>
#include<limits>
#include<boost/lexical_cast.hpp>
using std::clog; 
using std::endl;
namespace siesta{
	class ifstream;
	using namespace boost::units;

	static const divide_typeof_helper<atomic::length, si::femtosecond_unit>::type bohr_per_femtosecond;

	typedef divide_typeof_helper<nonsi::electron_volt_unit, nonsi::angstrom_unit>::type electronvolt_per_angstrom_unit;
	static const divide_typeof_helper<atomic::length, si::femtosecond_unit>::type electronvolt_per_angstrom;

	template<class Unit>
	struct vector : boost::array<double, 3>{
		//static vector from_value(boost::array<double, 3> const& arr){
		//	return 
		//}
		template<class Unit2>
		vector(vector<Unit2> const& v2){
			quantity<Unit> x(v2[0]*Unit2());
			quantity<Unit> y(v2[1]*Unit2());
			quantity<Unit> z(v2[2]*Unit2());
			(*this)[0] = x.value();
			(*this)[1] = y.value();
			(*this)[2] = z.value();
		}
		vector(
			quantity<Unit> const& x,
			quantity<Unit> const& y,
			quantity<Unit> const& z
		) : boost::array<double, 3>((boost::array<double, 3>){{x.value(), y.value(), z.value()}}){}
		/*
		friend vector<Unit> operator+(vector const& v1, vector const& v2){
			return vector<Unit>(
				(v1[0]+v2[0])*quantity<Unit>, 
				(v1[1]+v2[1])*quantity<Unit>, 
				(v1[2]+v2[2])*quantity<Unit>
			);
		}*/
		/*
		friend vector operator-(vector const& v1, vector const& v2){
			return vector(
				(v1[0]-v2[0])*quantity<Unit>, 
				(v1[1]-v2[1])*quantity<Unit>, 
				(v1[2]-v2[2])*quantity<Unit>
			);
		}*/
		vector operator-() const{
			return vector(
				-(*this)[0]*Unit(), 
				-(*this)[1]*Unit(), 
				-(*this)[2]*Unit()
			);
		}
		friend std::ostream& operator<<(std::ostream& os, vector const& v){
			return os << v[0] << " " << v[1] << " " << v[2] << " # cartesian vector with " << Unit() << " unit # ";
		}
	};
	template<class Unit1, class Unit2>
	typename multiply_typeof_helper<quantity<Unit1>, quantity<Unit2> >::type operator*(siesta::vector<Unit1> const& v1, siesta::vector<Unit2> const& v2){
		typedef typename multiply_typeof_helper<quantity<Unit1>, quantity<Unit2> >::type ret;
		return ret::from_value(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
	}
	template<class Unit1, class Unit2>
	typename siesta::vector<typename multiply_typeof_helper<Unit1, Unit2>::type> operator*(quantity<Unit1> const& q, siesta::vector<Unit2> const& v){
		typedef siesta::vector<typename multiply_typeof_helper<Unit1, Unit2>::type> ret;
		return ret(q*v[0]*Unit2(),q*v[1]*Unit2(),q*v[2]*Unit2());
	}
	template<class Unit1, class Unit2>
	typename siesta::vector<typename divide_typeof_helper<Unit2, Unit1>::type> operator/(siesta::vector<Unit2> const& v, quantity<Unit1> const& q){
		typedef siesta::vector<typename divide_typeof_helper<Unit2, Unit1>::type> ret;
		return ret(v[0]*Unit2()/q,v[1]*Unit2()/q,v[2]*Unit2()/q);
	}
	template<class Unit>
	vector<Unit> operator-(vector<Unit> const& v1, vector<Unit> const& v2){
		return vector<Unit>(
			(v1[0] - v2[0])*Unit(), 
			(v1[1] - v2[1])*Unit(), 
			(v1[2] - v2[2])*Unit()
		);
	}
	template<class Unit>
	vector<Unit> operator+(vector<Unit> const& v1, vector<Unit> const& v2){
		typedef siesta::vector<Unit> ret;
		return ret(
			(v1[0] + v2[0])*Unit(), 
			(v1[1] + v2[1])*Unit(), 
			(v1[2] + v2[2])*Unit()
		);
	}
	template<class Unit>
	siesta::vector<typename divide_typeof_helper<Unit, Unit>::type> versor(siesta::vector<Unit> const& v){
		return v/sqrt(v*v);
	}

	//typedef std::vector<std::pair<boost::array<double, 3>, unsigned> > coordinates;
	struct coordinates : public std::vector<std::pair<siesta::vector<nonsi::angstrom_unit>, unsigned> >{
		value_type const& operator[](size_type idx) const{return this->at(idx);}
		value_type& operator[](size_type idx){return this->at(idx);}
	};
	class string : boost::equality_comparable<string>{
		std::string impl_;
		void project_(){
			boost::to_lower(impl_);
			boost::erase_all(impl_, "_");
			boost::erase_all(impl_, ".");
		}
		public:
		string(std::string s) : impl_(s){project_();}
		string(char const c[]) : impl_(c){project_();}
//		string& operator=(std::string 
		operator std::string const&() const{return impl_;}
		bool operator==(string const& other) const{return (std::string const&)(*this)==(std::string const&)(other);}
		bool operator==(char const other[]) const{return operator==(siesta::string(std::string(other)));}
		friend std::ostream& operator<<(std::ostream& os, string const& self){
			return os<<self.impl_;
		}
		friend std::istream& operator>>(std::istream& is, string& self){
			std::string temp; is >> temp; self=temp; return is;
		}
	};
	class step{
		public:
		typedef std::vector<siesta::vector<electronvolt_per_angstrom_unit> > forces;
		long index;
		double E_KS_eV;
		double Etot_eV;
		coordinates outcoor;
		boost::array<double, 3> point_charge; //only present in modified version of siesta
		forces outforc;
		private:
		BOOST_STATIC_ASSERT( BOOST_VERSION > 104000 ); // / 100000 >= 1 && BOOST_VERSION / 100 % 1000 >= 40);
		std::ifstream& operator<<(std::ifstream& ifs){
			point_charge = (boost::array<double, 3>){{-1.,-1.,-1.}};
			using boost::spirit::qi::int_;
			using boost::spirit::qi::long_;
			using boost::spirit::qi::char_;
			using boost::spirit::qi::double_;
			using boost::spirit::qi::lexeme;
			using boost::spirit::qi::_1;
			using boost::spirit::qi::phrase_parse;

			using boost::spirit::ascii::space;
			using boost::phoenix::ref; 
			index = -1;
			E_KS_eV = std::numeric_limits<double>::quiet_NaN();
			Etot_eV = std::numeric_limits<double>::quiet_NaN();
			outforc.clear();
			outcoor.clear();
			//double Etot_eV;
			// seek begin of step
			//std::clog <<"looking for beggining of step" << std::endl;
			bool step_found=false;
			for(std::string line; getline((std::ifstream&)ifs, line, '\n'); ){ //std::clog << "line  " << line << std::endl;
				std::string::iterator first = line.begin(); // needed for boost 1.43
				step_found = phrase_parse(first, line.end(),
						("Begin elec-ion dyn. step = " >> long_[ref(index) = _1]),
						space);
				if(step_found){ //clog << "step " << idx << " found"<< endl;
					//ifs >> boost::ignore_line >> boost::ignore_line >> boost::ignore_line;
					break;
				}
			}
			if(not step_found){
				std::clog << "unexpected end of file" << std::endl;
				ifs.setstate(std::ios::failbit);
				return ifs;
			}
			for(std::string line; getline((std::ifstream&)ifs, line, '\n'); ){
				boost::algorithm::trim(line);
				std::string::iterator first = line.begin();
				if(line[0]=='o') std::clog <<"processing line '"<< line << "'" << std::endl;
				//if(line[0]=='s') std::clog <<"processing line "<< line << std::endl;
				if(line=="siesta: Atomic forces (eV/Ang):"){
					for(std::string line; getline((std::ifstream&)ifs, line, '\n'); ){
						double x, y, z; //int species; int idx;
						int idx;
						std::string::iterator first = line.begin();
						bool r = phrase_parse(first, line.end(),
								(
									int_[ref(idx)=_1] >> double_[ref(x) = _1] >> double_[ref(y)= _1] >> double_[ref(z)=_1] 
								),
								space);
						if(not r){
							//clog << "parsing failed on line " << line << endl;
							break;
						}
						outforc.push_back(forces::value_type(
							x*nonsi::electron_volt/nonsi::angstrom,
							y*nonsi::electron_volt/nonsi::angstrom,
							z*nonsi::electron_volt/nonsi::angstrom 
						));
					}
				}else if(line=="outcoor: Atomic coordinates (Ang):"){
					std::clog << "found coordinates" << std::endl;
					for(std::string line; getline((std::ifstream&)ifs, line, '\n'); ){
						double x, y, z; int species; int idx;
						std::string::iterator first = line.begin();
						bool r = phrase_parse(first, line.end(),
								(
									double_[ref(x) = _1] >> double_[ref(y)= _1] >> double_[ref(z)=_1] >> int_[ref(species)= _1] >> lexeme[+(char_ - ' ')] >> int_[ref(idx)=_1]
								),
								space);
						if(not r){
							break;
						}
						outcoor.push_back(coordinates::value_type(
							siesta::vector<nonsi::angstrom_unit>(
								x*nonsi::angstrom, 
								y*nonsi::angstrom,
								z*nonsi::angstrom
							), 
							species)
						);
					}
					std::clog << "found " << outcoor.size() << " atoms"<<std::endl;
				}else if(phrase_parse(first, line.end(), ("siesta: Etot    =" >> double_), space)){
					std::string::iterator first = line.begin();
					bool r = phrase_parse(first, line.end(), ("siesta: Etot    =" >> double_[ref(Etot_eV)=_1]), space);
					assert(r);
					std::clog << "found etot with " << Etot_eV << std::endl;
				}else if(phrase_parse(first, line.end(), ("siesta: E_KS(eV) =" >> double_), space)){
					std::string::iterator first = line.begin();
					bool r = phrase_parse(first, line.end(), ("siesta: E_KS(eV) =" >> double_[ref(E_KS_eV)=_1]),space);
					assert(r);
					std::clog << "found e_ks with " << E_KS_eV << std::endl;
				}else if(phrase_parse(first, line.end(), ("Point charge at" >> double_ >> double_ >> double_),space)){
					double x=-1, y=-1, z=-1;
					std::string::iterator first = line.begin();
					bool r = phrase_parse(first, line.end(), ("Point charge at" >> double_[ref(x)=_1] >> double_[ref(y)=_1] >> double_[ref(z)=_1]),space);
					assert(r);
					point_charge = (boost::array<double, 3>){{x,y,z}};
					//std::clog << "found point charge at " << x<<" " << y << " " << z <<  " while line was " << line << std::endl;
				}else if(line=="siesta:                 =============================="){
					break; //probably marks end of step;
				}
			}
			if(not outforc.empty() and not outcoor.empty() and (outforc.size() != outcoor.size())){
				std::clog << "warning: outforc.size() != outcoor.size() in siesta output, unexpected end of file" << std::endl;
			}
			return ifs;
		}
		friend class ifstream;
	};
	struct block{
	};
	struct ifstream;
	struct input : 
		std::map<std::string, std::string>
	{
		private: input(){}
		friend struct ifstream;
		public:
		std::map<std::string, std::string> block;
		
		// recognized keys
		double lattice_constant;
		// beg of recognized blocks
		boost::array< boost::array<double, 3>, 3> lattice_vectors;
		std::map<unsigned, std::pair<int, std::string> > chemical_species_label;
		coordinates atomic_coordinates_and_atomic_species;
		// end of recognized blocks

		void save(boost::filesystem::path p){
			boost::filesystem::ofstream ofs(p); ofs << std::setprecision(15);
			for(const_iterator it=begin(); it!=end(); ++it){
				ofs << it->first <<" " << it->second << std::endl;
			}
			for(std::map<std::string, std::string>::const_iterator it=block.begin(); it!=block.end(); ++it){
				ofs << "%block " << it->first << '\n' << it->second << "\n%endblock " << it->first << std::endl;
			}
			ofs << "%block lattice_vectors\n";
			for(unsigned i=0; i!=3; ++i){
				for(unsigned j=0; j!=3; ++j){
					ofs << lattice_vectors[i][j];
				}
				ofs << '\n';
			}
			ofs << "%endblock lattice_vectors\n";
			ofs << "%block chemical_species_label\n";
			for(std::map<unsigned, std::pair<int, std::string> >::const_iterator it=chemical_species_label.begin();
				it!=chemical_species_label.end(); ++it
			){
				ofs << it->first << "  " << it->second.first << "  " << it->second.second << "\n";
			}
			ofs << "%endblock chemical_species_label\n";
			ofs << "%block atomic_coordinates_and_atomic_species" << std::endl;
			for(unsigned i=0; i!=atomic_coordinates_and_atomic_species.size(); ++i){
				ofs
					<< atomic_coordinates_and_atomic_species[i].first[0] << ' '
					<< atomic_coordinates_and_atomic_species[i].first[1] << ' '
					<< atomic_coordinates_and_atomic_species[i].first[2] << ' '
					<< atomic_coordinates_and_atomic_species[i].second   << '\n'
				;
			}
			ofs << "%endblock atomic_coordinates_and_atomic_species" << std::endl;
		}
		void read(std::istream& is){
			for(std::string line; getline(is, line, '\n'); ){
				std::istringstream iss_line(line);
				//std::clog << "line is " << line << std::endl;
				std::string key; iss_line >> key;
				//std::clog << "key " << key<<" found" << std::endl;
				if(key[0]=='#'){
					continue; //comment;
				}else if(key=="%block"){
					std::string block_name;
					iss_line >> block_name; //block name;
					std::clog << "input block "<<key<<" found " <<std::endl;
					std::string block_value;
					for(std::string line; getline(is, line, '\n'); ){
						std::istringstream iss(line); std::string first; iss>>first;
						if(first=="%endblock"){
							//std::clog << "endblock found" << std::endl;;
							break;
						}
						block_value += line + '\n';
					}
					std::istringstream value_iss(block_value);
					if      (block_name=="chemical_species_label"){
						namespace qi = boost::spirit::qi;
						namespace px = boost::phoenix;
						namespace ascii = boost::spirit::ascii;
						unsigned idx; int atomic_number; std::string element_name;
						for(std::string line; getline(value_iss, line, '\n'); ){
							std::string::iterator first = line.begin(), last = line.end();
							if(
						        qi::phrase_parse(first, last,
					        		( qi::int_ >> qi::int_     >> qi::lexeme[+~qi::char_(' ')] ), ascii::space,
									idx,        atomic_number,  element_name
								)
							){
								chemical_species_label[idx] = std::pair<int, std::string>(atomic_number, element_name);
							}else{
								break;
							}
						}
					}else if(block_name=="lattice_vectors"){
						for(unsigned i=0; i!=3; ++i){
							for(unsigned j=0; j!=3; ++j){
								value_iss >> lattice_vectors[i][j];
							}
						}
					}else if(block_name=="atomic_coordinates_and_atomic_species"){
						using namespace boost::spirit::qi;
						using boost::phoenix::ref; 
						for(std::string line; getline(value_iss, line, '\n'); ){
							double x, y, z; int species; //int idx;
							std::string::iterator first = line.begin();
							bool r = phrase_parse(first, line.end(),
									(
										double_[ref(x) = _1] >> double_[ref(y)= _1] >> double_[ref(z)=_1] >> int_[ref(species)=_1]
									),
									space);
							if(not r){
								break;
							}
							atomic_coordinates_and_atomic_species.push_back(
								coordinates::value_type(
									siesta::vector<nonsi::angstrom_unit>(
										x*nonsi::angstrom,
										y*nonsi::angstrom,
										z*nonsi::angstrom
									), 
									species
								)
							);
						}
						std::clog << "found " << atomic_coordinates_and_atomic_species.size() << " atoms in input dump" << std::endl;
					}else{
						block[block_name]=block_value;
					}
				}else{
					std::string value; getline(iss_line, value, '\n');
					if(key == "lattice_constant"){
						lattice_constant = boost::lexical_cast<double>(value);
					}else{
						(*this)[key]=value;
					}
					//std::clog << "key " << key<<" found with value " << value << std::endl;
				}
			}
			// checks
			assert(atomic_coordinates_and_atomic_species.size()>0);
		}
		input(boost::filesystem::path const& p) try{
			boost::filesystem::ifstream ifs(p); if(not ifs) throw std::runtime_error("cannot open stream");
			assert(ifs);
			read(ifs);
		}catch(...){throw std::runtime_error("cannot load "+p.string());}
		input(std::istream& iss) try{
			read(iss);
		}catch(...){throw std::runtime_error("cannot read from stream");}
		input const& export_xyz(boost::filesystem::path const& p) const{
			boost::filesystem::ofstream ofs(p);
			ofs << this->atomic_coordinates_and_atomic_species.size() << "\n";
			ofs << "# comment line, volume ";
			{
				unsigned line = 0;
				for(
					coordinates::const_iterator it = atomic_coordinates_and_atomic_species.begin();
					it != atomic_coordinates_and_atomic_species.end(); ++it){
					ofs << "  " << chemical_species_label[it->second.second].second << " " << it->first[0] << " " << it->first[1] << " " << it->first[2];
					if(line == 0) {
						ofs
							    << " crystal_vector 1 " << s->atoms.cell().a(0)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom)
							    << " crystal_vector 2 " << s->atoms.cell().a(1)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom)
							    << " crystal_vector 3 " << s->atoms.cell().a(2)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom);
					}
				}
			}
			return *this;
		}
	};
	struct ifstream : std::ifstream{
		siesta::input input_dump;
		ifstream(boost::filesystem::path const& p) : std::ifstream(p.string().c_str()){
			if(not boost::filesystem::exists(p)) throw std::runtime_error("siesta output file " + p.string() + " can not be read");
			for(std::string line; std::getline(*this, line, '\n'); ){
				if(line=="************************** Dump of input data file ****************************"){
					std::string input;
					for(std::string line; std::getline(*this, line, '\n'); ){
						if(line=="************************** End of input data file *****************************"){
							std::clog << "end of input found" << std::endl;
							break;
						}
						input += line + '\n';
					}
					std::istringstream iss(input);
					input_dump = siesta::input(iss);
					break;
				}
			}
			assert(*this);
		}
		ifstream& operator>>(step& s){s<<(*this); return *this;}
		friend class step;
	};
}

#endif
#ifdef _SIESTA_HPP_TEST
int main(){
#if 0
	//test flexible_string;
	siesta::string key("AtomicCoordinatesAndAtomicSpecies");
	std::clog << key << std::endl;
	assert(key == "atomic_coordinates_and_atomic_species");
	key = "Hola";
	std::clog << key << std::endl;
	//test input read
	siesta::input i("siesta/test/siesta.fdf");
	std::cout << i.atomic_coordinates_and_atomic_species.size() << std::endl;
#endif
	//test short output
	siesta::ifstream si("siesta/test/short.sst"); // or siesta::output;
	si.input_dump.atomic_coordinates_and_atomic_species.size();
	for(siesta::step s; si>>s;){
		std::clog << "Etot " << s.Etot_eV << std::endl;
		std::clog << "step found" << std::endl;
	}
/*
	unsigned idx = 10;
	unsigned count = 0;
	std::ofstream ofs("siesta.dat");
	for(siesta::step s; si >> s;){
		unsigned n = s.outcoor.size();
		clog << "n atoms is " << n << std::endl;
		count++;
		ofs << count << " " << s.outforc[idx][1] << endl;
	}
*/
	return 0;
}
#endif
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:

