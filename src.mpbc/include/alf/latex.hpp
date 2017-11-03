#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ -I$HOME/prj -L$HOME/usr/lib -laspell $0.cpp -Wl,-rpath=$HOME/usr/lib -L$HOME/usr/lib `#-std=c++0x` -Wall -Wextra `#-Wfatal-errors` -Wno-unused-variable -Wno-unused-parameter -Wno-ignored-qualifiers -lhunspell -laspell -I$HOME/usr/include   -lboost_regex -lboost_system -lboost_filesystem  -D_TEST_LATEX_HPP -o ./$0.x && ./$0.x $@ && (acroread latex.pdf&)
rm -rf $0.cpp $0.x
exit
#endif
#ifndef LATEX_HPP
#define LATEX_HPP
#include<string>
#include<sstream>
#include<vector>
#include<boost/filesystem/fstream.hpp>
#include<boost/filesystem/convenience.hpp>
#include<boost/format.hpp>
#include<map>
#include<set>
#include<complex>
#include<iostream> //clog
#include<boost/lexical_cast.hpp>
#include<boost/tuple/tuple.hpp>
#include<boost/regex.hpp>
#include"boost/named.hpp"
#include<iomanip>
#include<boost/math/special_functions/fpclassify.hpp> //(boost::math::isinf)
#include <boost/algorithm/string.hpp> // boost::replace_last
#include<map>
#include<boost/numeric/interval.hpp>
#include<boost/multi_array.hpp>
#include<boost/process.hpp>
#include<boost/algorithm/string/classification.hpp>
#include<boost/algorithm/string/split.hpp>
#include<boost/logic/tribool.hpp>
#include<boost/tokenizer.hpp>
#include<boost/units/systems/si.hpp>
//#include "alf/hunspell.hpp" // needs -lhunspell
#include "alf/aspell.hpp" // needs -laspell
using boost::lexical_cast;

#ifndef LATEX_UNITS_HPP
#define LATEX_UNITS_HPP
#include<boost/units/quantity.hpp>
#include<boost/units/physical_dimensions.hpp>
namespace boost{
namespace units{
using std::string;
template<class Dimension> string name                 (){return "";}

template<> string name<energy_dimension              >(){return "Energy";}
template<> string name<force_dimension               >(){return "Force";}
template<> string name<velocity_dimension            >(){return "Velocity";}
template<> string name<volume_dimension              >(){return "Volume";}
template<> string name<length_dimension              >(){return "Length";}
template<> string name<pressure_dimension            >(){return "Pressure";}
template<> string name<temperature_dimension         >(){return "Temperature";}
template<> string name<heat_capacity_dimension       >(){return "Entropy";}
template<> string name<time_dimension                >(){return "Time";}
template<> string name<mass_density_dimension        >(){return "Mass Density";}
template<> string name<thermal_conductivity_dimension>(){return "Thermal Conductivity";}
template<> string name<wavenumber_dimension          >(){return "Wavenumber";}
template<> string name<
	derived_dimension<
		length_base_dimension,  -3,
		time_base_dimension,     1,
		current_base_dimension,  1
	>::type
>(){return "Charge Density";}

}}

namespace boost{
namespace units{
using std::string;
inline string latexify(string const& units_s){
	using namespace boost;
	std::string ret = units_s;
	ret=regex_replace(ret,
		regex("\\^\\(*([+,-]*\\d{1,}/*\\d{0,})\\)*", boost::regex::perl),
		"^{$1}"
	);
	ret=regex_replace(ret,
		regex("(\\s)", regex::perl),
		"\\\\, "
	);
	boost::regex subscripts("_(\\w+)", boost::regex::icase|boost::regex::perl);	
	ret=regex_replace(ret,subscripts,"_{\\\\rm $1}");
	return /*"\\mathrm{"+*/ ret /*+"}"*/
	;
}

/// Transforms usual Boost unit symbol_string into usual LaTeX output
template<class Dimension, class System>
inline std::string latex_string(unit<Dimension, System> const& r){ // , boost::disable_if< boost::is_same<System, si::system>, void >* dummy = 0){
	return latexify(symbol_string(r));
}

//> this is now in the namespace 'si' below
//template<class Dimension>
//inline std::string latex_string(unit<Dimension, si::system> const& r){
//	return "\\rm{" + latexify(symbol_string(r)) + "}";
//}

template<class System>
inline std::string latex_string(unit<dimensionless_type, System> const&){
	return "";
}


}}

#include<boost/units/systems/si.hpp>
#include"boost/units/systems/atomic.hpp"
namespace boost{
namespace units{
namespace si{
//template<class Unit>
std::string latex_string(si::femto_scaled<si::time>::unit const&){
	return "\\ensuremath{\\rm{fs}}"; //+latex_string(si::time());
}
}}}
namespace boost{
namespace units{

using std::string;
template<> string name<atomic::wavefunction_amplitude::dimension_type>(){return "Wavefunction Amplitude";}

namespace si{
template<class Dimension>
inline std::string latex_string(unit<Dimension, si::system> const& r){
	return "\\rm{"+latexify(symbol_string(r))+"}";
}
inline std::string latex_string(unit<dimensionless_type, si::system> const&){
	return "";
}
}}}

#include<boost/units/systems/cgs.hpp>
namespace boost{
namespace units{
namespace cgs{
template<class Dimension>
inline std::string latex_string(unit<Dimension, cgs::system> const& r){
	return "\\rm{"+latexify(symbol_string(r))+"}";
}
}}}


namespace boost{
namespace units{
namespace nonsi{
std::string latex_string(nonsi::electron_volt_unit const&){
	return "\\ensuremath{\\rm{eV}}";
}
}}}
#endif

/// Automatic generation of LaTeX code for output
namespace latex{
	using std::string;
	using std::endl;
	using std::clog;
	using namespace boost::filesystem;
	using namespace boost::units;
	struct ostream;
	struct language : virtual std::string{
		language(){}
		language(std::string const& s) : std::string(s){}
		virtual void operator>>(ostream& os) const;
	};
	struct english : language{
		english(std::string const& s) : std::string(s){}
		virtual void operator>>(ostream& os) const;
	};
	struct error : std::runtime_error{
		error(std::string const& what) : std::runtime_error(what){}
	};
	struct undefined_control_sequence : error{
		undefined_control_sequence(std::string const& what) : error(what){}
	};
	class document{};
	class article : document{};
	static std::map<void*, std::string> name;
	struct environment;
	/// LaTeX (or indirect PDF) output stream, it can feed from specific LaTeX object of from strings with LaTeX code
	struct ostream : protected std::ostringstream{
		static std::map<string, string> unicode_translation;
		protected:
		using std::ostringstream::str;
		public:
		std::set<path> temp_files_;
		std::set<string> packages_;
		std::map<string, string> fields_;
		string title;
		string author;
		std::ostringstream preamble_;
		typedef std::ostringstream base;
		path p_;
		ostream(){}
		ostream(path const& p) : p_(p) {};

		friend ostream& operator<<(ostream& os, bool   const& b){
			if(b==true){
				((base&)(os))<< "\\ding{51}"; // pifont::ding{51} // amsmath::checked unicode::2713 http://en.wikipedia.org/wiki/List_of_Unicode_characters#Dingbats
			}else{
				((base&)(os))<< "\\ding{55}"; // ding: pifont package unicode::2717 http://en.wikipedia.org/wiki/List_of_Unicode_characters#Dingbats
			}
			 return os;
		}
		
		friend ostream& operator<<(ostream& os, boost::logic::tribool const& tb){
			if(tb){
				((base&)(os))<< "\\mvchr{86}"; // marvosym::Checkedbox pifont::CheckedBox  unicode::2611 ballot box with check ☑
			}else if(!tb){
				((base&)(os))<< "\\mvchr{88}"; // marvosym::Crossedbox pifont::XBox unicode::2612 ballot box with x ☒
			}else{
				((base&)(os))<< "\\mvchr{79}"; //"\\Square"; //  pifont::Square unicode: 2610 ballot box ☐
			}
			return os;
		}
		friend ostream& operator<<(ostream& os, int      const& s){((base&)(os))<< s; return os;}
		friend ostream& operator<<(ostream& os, unsigned const& s){((base&)(os))<< s; return os;}
		friend ostream& operator<<(ostream& os, double const& d); 
		friend ostream& operator<<(ostream& os, string const& s){
			((base&)(os))<< s; return os;
		}
		friend ostream& operator<<(ostream& os, language const& l){ l >> os; return os; };
		friend ostream& operator<<(ostream& os, const char* const c){
			//string c_str(c);
			//if(c_str.size()>80) 
			//	os << english(c_str);
			//else
			((base&)(os)) << c; 
			return os;
		}
		friend ostream& operator<<(ostream& os, boost::filesystem::path const& p){os<<"\\url{file:" << p.string() << "}";return os;}
		friend ostream& operator<<(ostream& os, environment const& env);
		template<class Unit>
		friend ostream& operator<<(ostream& os, boost::units::quantity<Unit> const& q);
		~ostream(){
			flush(*this);
			for(std::set<path>::iterator it=temp_files_.begin(); it!=temp_files_.end(); ++it){
				clog << "please delete 'rm " << it->string() << "'" << std::endl;
				//boost::filesystem::remove(*it);
			}
		}
		void declare_unicode_character(string const& code, string const& symbol){
			preamble() << ("\\DeclareUnicodeCharacter{\""+ code +"}{"+ symbol +"}\n");
		}
		static ostream& flush(ostream& os){
			std::clog << "generating " << os.p_.string() << std::endl;
			//std::string test_line = "! LaTeX Error: File `aaa.tex' not found.";
			//boost::regex file_not_found("LaTeX Error File `(.*)' not found");//("! LaTeX Error: File `(.*)' not found.");
			
			(base&)os<<std::flush;
			if(os.p_==path("")) return os;
			path p_tex;
			path p_pdf;
			if(os.p_.extension()==".tex"){
				p_tex = os.p_;
			}else if(os.p_.extension()==".pdf"){ // note: latex auxiliary files can't be invisible
				p_tex = os.p_; p_tex.replace_extension(".tex");
				//p_tex = string(os.p_.stem()) + string(".tex"); //ok for boost 1.42, not ok boost 1.46
				//p_tex = os.p_.stem().string() + string(".tex"); //ok for boost 1.46, not ok for boost 1.42
			}else{ //< possibly add dvi and ps, dvi to text also is interesting
				throw std::runtime_error("do not know what to do with given extension for LaTeX output"); //"+os.p_.extension() +" .string() 
			}
			{
				boost::filesystem::ofstream ofs(p_tex);
				ofs << "\\documentclass[12pt,fleqn]{tufte-handout}\n"; //{article}\n";
				ofs << "\\usepackage["
					<< "compatibility=off, "
					<< "font={small}]{caption}\n"; //for caption* (not numbered), compatibility=off to ensure this
				//ofs <<
					//"\\usepackage["
					//"mathletters"// autogenerated" //, postscript"
					//"]{ucs}\n"; // to do something useful with unicode characters in document/cpp source
				ofs << "\\usepackage{comment}";
				ofs << "\\usepackage["
					<< "utf8x"
					<< "]{inputenc}\n";
				ofs << "\\usepackage{mathrsfs}\n"; // for mathscr
				//ofs << "\\usepackage{cool}\n";
				ofs << "\\usepackage{listings}\n";
				ofs << "\\usepackage{dcolumn}\n";
				ofs << "\\usepackage{xcolor}\n";   //included in tufte-latex but useful to convert document to plain article
				ofs << "\\usepackage{geometry}\n"; //idem 
				ofs << 
					"\\definecolor{darkgray}{rgb}{0.95,0.95,0.95}\n "
					"\\lstset{language={[GNU]C++}}\n "
					"\\lstset{backgroundcolor=\\color{darkgray}}\n"
					"\\lstset{numbers=left, numberstyle=\\tiny, stepnumber=2, numbersep=5pt}\n"
					//"\\lstset{keywordstyle=\\color{red}\\bfseries\\emph}\n"
				;
				ofs <<"\\hypersetup{pdfencoding=unicode}\n"; //hyperref options set in w/hypersetup because tufte loads hyperref
				// ofs << "\\usepackage[pdfborder={0 0 0}, pagebackref=true]{hyperref}\n"; //not compatible with tufte-*
				// ofs<<"\\usepackage{courier}\n"; //look like ms courier 
				ofs << "\\newcommand{\\unicodesub}[1]{_{#1}}\n"; // ₀
				ofs << "\\newcommand{\\unicodewide}[1]{#1}\n";
				//ofs << "\\newcommand{\\unicodesuper}[1]{^{#1}}\n";
				os.declare_unicode_character("00B2", "^2");                                       // ²
				ofs << "\\DeclareUnicodeCharacter{\"BD}{\\ensuremath{\\text{\\textonehalf}}}\n";   // ½
				ofs << "\\DeclareUnicodeCharacter{\"035D}{\\ensuremath{\\phi}} \n";                // �  math phi
				os.declare_unicode_character("03B1", "\\alpha");                                   // α theta, sometimes translated as \texttheta, result of copying from pdf avoids ! LaTeX Error: Command \textalpha unavailable in encoding T1.
				os.declare_unicode_character("03B2", "\\beta");                                    // β theta, sometimes translated as \texttheta, result of copying from pdf avoids ! LaTeX Error: Command \textalpha unavailable in encoding T1.		
				os.declare_unicode_character("03B8", "\\theta");                                   // θ theta, sometimes translated as \texttheta, result of copying from pdf avoids ! LaTeX Error: Command \textalpha unavailable in encoding T1.
				ofs << "\\DeclareUnicodeCharacter{\"212B}{\\mathrm{\\mathring{A}}} \n";            // Å \\ensuremath{\\rm{\\AA}}
				ofs << "\\DeclareUnicodeCharacter{\"2329}{\\ensuremath{\\langle}}\n";              // 〈
				ofs << "\\DeclareUnicodeCharacter{\"232A}{\\ensuremath{\\rangle}}\n";              // 〉
				ofs << "\\DeclareUnicodeCharacter{\"1D62}{\\ensuremath{_i}}\n";                    // ⱼ 
				ofs << "\\DeclareUnicodeCharacter{\"2044}{\\over}\n";                              // ⁄ fraction slash
				ofs << "\\DeclareUnicodeCharacter{\"2C7C}{\\ensuremath{_j}}\n";                    // ⱼ 
				ofs << "\\DeclareUnicodeCharacter{\"2148}{\\unichar{105}}\n";                      //
				ofs << "\\DeclareUnicodeCharacter{\"1D4E9}{\\ensuremath{\\mathcal{Z}}}\n";         // ⱼ 
				ofs << "\\DeclareUnicodeCharacter{\"1D4D4}{\\ensuremath{\\mathcal{E}}}\n";         // ⱼ 
				ofs << "\\DeclareUnicodeCharacter{\"1D4D4}{\\ensuremath{\\mathcal{E}}}\n";         // ⱼ 
				os.declare_unicode_character("1D609", "_\\text{B}");                                // �
				os.declare_unicode_character("1D624", "_\\text{c}");                               // 𝘤
				os.declare_unicode_character("1D6FA", "\\Omega");                                  // 𝛺
				os.declare_unicode_character("1D6FD", "\\beta");                                  // 𝛽
				os.declare_unicode_character("1D707", "\\mu");                                     // �
				os.declare_unicode_character("1D70E", "\\sigma");
				os.declare_unicode_character("1D708", "\\nu");                               // �
				ofs << "\\DeclareUnicodeCharacter{\"2192}{\\ensuremath{\\rightarrow}}\n";          // → rightwards arrow
				ofs << "\\DeclareUnicodeCharacter{\"2799}{\\ensuremath{\\to}}\n";                  // ➙ rightwards arrow
				ofs << "\\DeclareUnicodeCharacter{\"27F6}{\\ensuremath{\\longrightarrow}}\n";      // ⱼ 
				ofs << "\\DeclareUnicodeCharacter{\"2423}{\\ensuremath{\\quad}}\n";                // ␣ open box graphic for space
				ofs << "\\DeclareUnicodeCharacter{\"FF08}{\\bigg(}\n";                             // （ open wide parethesis
				ofs << "\\DeclareUnicodeCharacter{\"FF09}{\\bigg)}\n";                             // （ close wide parethesis
				ofs << "\\DeclareUnicodeCharacter{\"FF5B}{\\left\\{}\n";                           // ｛ open wide curly bracket
				ofs << "\\DeclareUnicodeCharacter{\"FF5D}{\\right\\}}\n";                          // ｝ close wide curly bracket
				ofs << "\\DeclareUnicodeCharacter{\"FFFD}{\\ensuremath{\\mathcal{E}}}\n";          //
				ofs << "\\DeclareUnicodeCharacter{\"1D700}{\\ensuremath{\\varepsilon}}\n";         //
				ofs << "\\DeclareUnicodeCharacter{\"1D714}{\\ensuremath{\\omega}}\n";              //
				ofs << "\\DeclareUnicodeCharacter{\"1D70B}{\\ensuremath{\\pi}}\n";                 //
				ofs << "\\usepackage{amsmath}\n"; //for begin{cases}
				ofs << "\\usepackage{marvosym}\n"; //for checkedbox etc
				ofs << "\\usepackage{wasysym}\n";
				ofs << "\\usepackage{pifont}\n";
				ofs << "\\usepackage{bbold}\n"; //bbold for mathbb
				ofs << "\\usepackage{textcomp}\n"; //for capitalring
				ofs << "\\usepackage{color}\n";
				//ofs << "\\usepackage[left=1cm, right=1cm, top=1.7cm, bottom=1.7cm]{geometry}\n"; //not compatible with tufte-*
				//ofs << "\\geometry{left=1.5cm, right=10cm, marginparwidth = 7cm}\n"; //right=1cm, top=1.7cm, bottom=1.7cm
				ofs << "\\geometry{left=1.1cm, right=8.5cm, marginparwidth = 6.5cm}\n"; //right=1cm, top=1.7cm, bottom=1.7cm
				// include pdfcomment  (pdf annotations), savesymbol used to make it compatible with tufte-latex!
				ofs << 
					"\\usepackage{savesym}\n"
					"\\savesymbol{marginnote}\n"
					"\\usepackage{pdfcomment}\n"
					"\\restoresymbol{pdfcomm}{marginnote}\n"
					"\\savesymbol{Cross}\n" //avoid name conflict between cool/bbm
					"\\usepackage{cool}\n" // functions with natural names, e.g. Erf
				;
				//ofs << "\\usepackage{pdfcomment}\n"; //not compatible with tufte-* important
				ofs << "\\usepackage[]{attachfile2}\n";
				ofs << "\\usepackage[]{embedfile}\n";
				ofs << "\\usepackage[]{pdfmarginpar}\n";
				ofs << "\\pdfstringdefDisableCommands{%\n";
				ofs << "     \\let\\n\\textLF\n";
				ofs << " }\n"; //to allow linebreaks in pdf-comments
//				ofs << "\\long\\def\\ASYaligned(#1,#2)(#3,#4)#5#6#7{\\leavevmode% \n"
//					" \\setbox\ASYbox=\\hbox{#7}% \n"
//					" \\setbox\\ASYbox\\hbox{\\ASYdimen=\\ht\\ASYbox%\n "
//					" \advance\ASYdimen by\dp\ASYbox\kern#3\wd\ASYbox\raise#4\ASYdimen\box\ASYbox}% \put(#1,#2){#5\wd\ASYbox 0pt\dp\ASYbox 0pt\ht\ASYbox 0pt\box\ASYbox#6}}% \long\def\ASYalignT(#1,#2)(#3,#4)#5#6{% \ASYaligned(#1,#2)(#3,#4){% \special{pdf:q #5 0 0 cm}% }{% \special{pdf:Q}% }{#6}} \long\def\ASYalign(#1,#2)(#3,#4)#5{\ASYaligned(#1,#2)(#3,#4){}{}{#5}} \def\ASYraw#1{#1} \n";
				for(std::set<string>::const_iterator it=os.packages_.begin(); it!=os.packages_.end(); ++it){
					ofs<<"\\usepackage{"<<*it<<"}\n";
				}
				ofs<<os.preamble().str()<<"\n";
				//ofs << "\\usepackage[active,tightpage]{preview}\n"
				//	"\\PreviewEnvironment{tikzpicture}\n";
				//ofs<<"\\usepackage{cmbright}\n";
				//ofs<<"\\usepackage[math]{iwona}\n";
				//ofs<<"\\usepackage[math]{anttor}\n";
				//ofs<<"\\usepackage{pxfonts}\n"; //ofs<<"\\usepackage{mathpazo}\n"; //ofs<<"\\usepackage{mathpple}\n";
				//ofs<<"\\usepackage[varg]{txfonts}\n"; //\usepackage{mathtime} //\usepackage{mathptmx} //\usepackage{mbtimes} //times like
				//ofs<<"\\usepackage{arev}\n"; // blocky
				//ofs<<"\\usepackage[charter]{mathdesign}\n"; //very nice
				//ofs<<"\\usepackage{fourier}\n"; //nice
				//ofs<<"\\renewcommand{\\rmdefault}{cmss}\n";
				//ofs<<"\\usepackage[firstpage]{draftwatermark}\n\\SetWatermarkLightness{ 0.95 }\n\\SetWatermarkScale{ 5 }\n"; //can't use this because it doesn't work in single page documents
				string watermark ; //= "DRAFT";
				ofs<<
					"\\usepackage{graphicx}\n"
					"\\usepackage{type1cm}\n"
					"\\usepackage{eso-pic}\n" //for watermark, other packages fail to put a water mark in the last page
					"\\usepackage{color}\n"
					"\\makeatletter\n"
					"\\AddToShipoutPicture{%\n"
					"   \\setlength{\\@tempdimb}{.5\\paperwidth}%\n"
					"   \\setlength{\\@tempdimc}{.5\\paperheight}%\n"
					"   \\setlength{\\unitlength}{1pt}%\n"
					"   \\put(\\strip@pt\\@tempdimb,\\strip@pt\\@tempdimc){%\n"
					"	\\makebox(0,0){\\rotatebox{45}{\\textcolor[gray]{0.95}%\n"
					"	{\\fontsize{6cm}{6cm}\\selectfont{" + watermark + "}}}}%\n"
					"}%\n"
					"}\n"
					"\\makeatother";
				//ofs << "\\usepackage{natbib}\n";
				//ofs << "\\bibliographystyle{plainnat}\n";
				ofs << "\\usepackage{backref}\n";
				ofs << "\\usepackage{lastpage}\n";
				ofs << "\\usepackage{fancyhdr}\n\\setlength{\\headheight}{15.2pt}\n\\pagestyle{fancy}\n";
				ofs << "\\usepackage{datetime} \\ddmmyyyydate\n"; //{\\color{gray}{\\currenttime}}~ \\usdate
				//ofs << "\\lfoot{Alfredo Correa \\href{mailto:correaa@llnl.gov}{\\nolinkurl{<correaa@llnl.gov>}}}\n\\cfoot{\\thepage \\color{gray}{/\\pageref{LastPage}}}\n\\rfoot{\\today~\\color{gray}{\\currenttime}}\n"; //\\color{gray}
				//ofs << "\\renewcommand{\\headrulewidth}{0.0pt}\n\\renewcommand{\\footrulewidth}{0.4pt}\n";
				ofs << "%\\usepgfplotslibrary{external}\n"
					  "%\\tikzexternalize% activate externalization!\n";
				ofs << "\\widowpenalty=10000\n"; //up to 10000 from http://en.wikibooks.org/wiki/LaTeX/Page_Layout#Widows_and_orphans
				ofs << "\\clubpenalty=10000\n"; //up to 10000  or from http://xpt.sourceforge.net/techdocs/language/latex/latex35-OrphanSpaceControl/ar01s04.html
				ofs << "\\setlength{\\parskip}{3ex plus 2ex minus 2ex}\n"; //rubber band between paragraphs
				if(os.title!="") ofs<<"\\title{"<<os.title<<"}\n";
				if(os.title!="") ofs<<"\\author{"<<os.author<<"}\n";
				ofs <<
					"\\begin{document}\n"
						<< os.str() << "\n\n"
					"\\end{document}\n"
				;
			}
			if(os.p_.extension()==".pdf"){
				clog << "generating "<< p_tex.string() << " ... "<< std::endl;
				std::string exec = boost::process::find_executable_in_path("pdflatex"); //"pdflatex"; 
				std::vector<std::string> args; 
				args.push_back(exec);
				args.push_back("-shell-escape"); 
				//args.push_back("-interaction");
				//args.push_back("nonstopmode"); //("batchmode");
				//args.push_back("-halt-on-error");
				args.push_back(p_tex.string());
				boost::process::context ctx;
				ctx.environment = boost::process::self::get_environment();
				ctx.stdin_behavior = boost::process::capture_stream();
				ctx.stdout_behavior = boost::process::capture_stream();//boost::process::silence_stream();
				for(unsigned i=0; i<2; ++i){
					boost::process::child c = boost::process::launch(exec, args, ctx); 
					boost::process::pistream &is = c.get_stdout();
					boost::process::postream &os = c.get_stdin();
					std::string line;
					std::string error_type;
					std::string error_context;
					while (std::getline(is, line)){
						if(line[0]=='!'){   
							error_type=line;
							clog << "warning: LaTeX error " << error_type <<  " found" << std::endl;
							{
								//std::string test_line = "! LaTeX Error: File `aaa.tex' not found.";
								//boost::regex test_rx("! LaTeX Error: File `aaa.tex' not found.");
								//std::clog << boost::regex_replace(test_line, test_rx, "$1");
							}
							if(error_type=="! Undefined control sequence."){
								std::getline(is, error_context);
								std::vector<std::string> v;
								boost::algorithm::split(v, error_context, boost::algorithm::is_any_of(" "));
								std::clog << "error in " << p_tex.string() <<" : undefined " << v[v.size()-2] << std::endl;
								os  
									<< "I{\\color{red}\\verb+" 
									<< v[v.size()-2]
									<< "+}" << std::endl
								;
							//}else if(boost::regex_match(error_type, file_not_found)){
							//	{std::ofstream ofs("missing.tex"); ofs << "{\\color{red}\\verb+<missing file>+}" << std::endl;}
							//	clog << "warning: file not found "<< error_type << std::endl; 
							//	os << "missing.tex" << std::endl;
							}else{
								std::string error_inserted;
								std::getline(is, error_inserted);
								std::clog << "error in "<< p_tex.string() <<" : inserted " << error_inserted << std::endl;
								os << " " << std::endl //ignore error
								;
							}
						}
					}
				}//for
			}
			std::clog << "generated " << os.p_.string() << std::endl;
			return os;
		}
		ostream& usepackage(string const& s){
			packages_.insert(s); 
			return *this;
		}
		static ostream& endl(ostream& os){
			//os<<"\n\n"; 
			//return flush(os);
			return flush(os<<"\n\n");
		}
		std::ostringstream& preamble(){return preamble_;}
		string preamble() const{return preamble_.str();}
		typedef std::ostream&(*std_manipulator)(std::ostream&);
		ostream& operator<<(std_manipulator manip){
			if(manip==(std::ostream& (*)(std::ostream&))std::endl){ // need cast due to some overload
				return ostream::endl(*this); //this can avoid lots of confusion
			}
			if(manip==(std::ostream& (*)(std::ostream&))std::flush){
				return ostream::flush(*this);
			}
			return *this;
		}
	};
	void language::operator>>(ostream& os) const{
		os << (std::string const&)(*this);
	}
	void english::operator>>(ostream& os) const{
		std::string corrected(*this);
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
		tokenizer tok(*this, boost::char_separator<char>(" .,;:[]~!@#$%^&*()_+=-`\"\'?<>/"));
		//using aspell::speller;
		aspell::speller c;
		for(tokenizer::iterator it = tok.begin(); it != tok.end(); ++it){
		    std::cout << *it <<" "; 
			if(it->size()>2 and (*it)[0]!='\\' and boost::regex_match(*it, boost::regex("^[a-z]*$") ) ){
				if(c(*it)){
				}else{
					std::string comment; 
					std::list<string> suggestions = c.suggest(*it);
					for(std::list<string>::const_iterator it_suggest= suggestions.begin(); it_suggest!=suggestions.end(); ++it_suggest){
						comment += (*it_suggest) + string(boost::next(it_suggest)!=suggestions.end()?"\\n ":""); //string((i!=suggestions.size()-1)?"\\n ":"");
					}
					boost::algorithm::replace_all(corrected, 
						" "+(*it)+" ", 
						" \\pdfmarkupcomment[author={not \"" + *it + "\" in dictionary}, subject={English dictionary}, color=Red, opacity=0.5, markup=Squiggly]{" + *it + "}{" + comment +"}\n " 
					);
				}
			}
			std::cout<< std::endl;
		}
		os << corrected.c_str();
	}
	ostream& endl(ostream& os){
		return ostream::endl(os);
	}
	ostream& flush(ostream& os){
		return ostream::flush(os);
	}
	struct par{
		friend ostream& operator<<(ostream& os, par const&){
			return os<<"\n\n\\par ";
		}
	} par;
	struct noindent{
		friend ostream& operator<<(ostream& os, noindent const&){
			return os<<"\\noindent \n";
		}
	} noindent;
	struct newpage{
		friend ostream& operator<<(ostream& os, newpage const&){
			return os<<"\\newpage\n";
		}
	} newpage;
	template<class Command>
	struct nullary{
//		friend ostream& operator<<(ostream& os, Command const&); //{return os<<"\\"<<Command::name()<<"\n";}
	};
	template<class Command>
	struct unary{
		string code_;
		string postscript_;
		unary(string code) : code_(code){}
		friend ostream& operator<<(ostream& os, Command const& c){
			return os<<"\\"<<Command::name()<<"{"<<c.code_<<"}" << c.postscript_;
		}
	};
	struct maketitle : nullary<maketitle>{
		static string name(){return "maketitle";}
	} maketitle;
	struct bibliography : unary<bibliography>{
		bibliography(string bibfile) : unary<bibliography>(bibfile){
			//assert(boost::filesystem::exists(bibfile+".bib")); //need to link to boostfileystem
			
			postscript_ = "\\attachfile[appearance=false, print=false, subject=Bibliography, description=BibTeX format, mimetype=text/plain]{" + bibfile + ".bib}";
		}
		static string name(){return "bibliography";}
	};
	struct bibliographystyle : unary<bibliographystyle>{ 
		bibliographystyle(string style /*eg. unsrtnat*/) : unary<bibliographystyle>(style){}
		static string name(){return "bibliographystyle";}
	};
	struct ref : unary<ref>{ // ref def
		ref(std::string const& s) : unary<ref>(s){}
		ref(environment const* ptr) : unary<ref>(boost::lexical_cast<string>(ptr)){}
		static string name(){return "ref";}
	};
	template<class Command>
	ostream& operator<<(ostream& os, nullary<Command> const&){return os<<"~\\"<<Command::name()<<"\n";} // ~ to avoid annoying error: ! LaTeX Error: There's no line here to end.
	ostream& operator<<(ostream& os, nullary<class maketitle> const&){
		assert(os.title!="");
		return os<<"\\"<<maketitle::name()<<"\n"; 
	}
	struct title{
		std::string title_;
		std::string author_;
		title(
			std::string const& t, 
			std::string const& a
		) : title_(t), author_(a){}
		friend ostream& operator<<(ostream& os, title const& t){
			os.title = t.title_;
			os.author = t.author_;
			os << maketitle;
			return os;
		}
	};
	struct newline : nullary<newline>{
		static string name(){return "newline";} 
	} newline;
	struct math{
		string code_;
		template<class T>
		math(T* tp) : code_(latex::name[(void*)tp]){}
		math(const char* code)   : code_(code){if(code_==" ") code_="\\quad";}
		math(char const code)    : code_(boost::lexical_cast<string>(code)){if(code_==" ") code_="\\:";} //other options: http://chenfuture.wordpress.com/2008/03/22/math-spacing-and-length-units/
		math(string const& code) : code_(code){}
		math(int code) : code_(boost::lexical_cast<string>(code)){}
		template<class Unit>
		math(boost::units::quantity<Unit> const& q) : code_("0"){
			if(q.value()==0){return;}
			//std::ostringstream oss; oss << boost::format("%1$1.3d") % q.value();
			code_ = 
				"{"
				"\\href{http://www.wolframalpha.com/input/?i=" + to_string(q.value())/*oss.str()*/ + "*(" + name_string(Unit()) + ")}\n"
				+ "{" + 
				math(q.value()).code_ + "~" /*~\\cdot*/ + latex_string(Unit())
				+"}"
				"}"
			;
		}
		math(double const& d){
			std::ostringstream oss;
			oss << /*boost::format("%1$9.4d") %*/ d ; // << "\\ddots"; //"{\\color{gray}{'}}";
			code_ = oss.str();
			using namespace boost;
			if(not find_last(code_,".")){
				if(find_last(code_, "e")){
					code_ = regex_replace(
						code_,
						regex(
							"e", 
							boost::regex::icase|boost::regex::perl
						), 
						".e"
					);
				}else{
					code_ += ".";
				}
			}
			code_ = regex_replace(
				code_,
				regex(
					"e\\+*(-*)0*(\\d{1,})", 
					boost::regex::icase|boost::regex::perl
				), 
				//"{\\\\color{gray}'}"
				"{\\\\scriptstyle \\\\times {10}^{$1$2}}"
			);
			code_ = regex_replace(
				code_,
				regex(
					"\\+0{1,}(\\d{1,})", 
					boost::regex::icase|boost::regex::perl
				), 
				"$1"
			);
		}
		math(std::complex<double> const& c){
			{
				std::ostringstream oss; oss << boost::format("%1$1.3d") % c.real(); oss.str();
				code_ = math(c.real()).code_;
			}{
				if(c.imag()==0){
					code_ = "(" + code_ + "_{ \\color{gray}{+ \\mathrm{i} 0.}})";
				}else{
					// std::ostringstream oss; oss << boost::format("%1$1.3d") % c.imag();
					code_ = "(" + code_ + ((c.imag()>0.)?" + ":"") + math(c.imag()).code_ + "\\mathrm{i})";  //oss.str();
				}
			}
		}
		friend ostream& operator<<(ostream& os, math const& self){return os<<self.code();}
		math operator()(math const& other) const{
			return math(this->code_ + "("+ other.code_ +")");
		}
		math operator()(math const& other1, math const& other2) const{
			return math(this->code_ + "("+ other1.code_ +", "+ other2.code_ + ")");
		}
		math operator=(math const& other) const{ //maybe has to be operator==
			return math( this->code_ + " = " + other.code_);
		}
		math operator|(math const& other) const{
			return math( "\\left." + this->code_ + "\\right|_{" + other.code_ + "}");
		}
		math operator,(math const& other) const{
			return math( this->code_ + " , " + other.code_);
		}
		math operator[](string const& subname){
			return math( "{"+this->code_ + "}_{\\rm{"+subname+"}}");
		}
		math operator[](int subindex){
			return math("{"+this->code_+"}_{"+lexical_cast<string>(subindex)+"}");
		}
		math operator[](math const& other) const{
			return math("{"+this->code_+"}_{"+other.code_+"}");
		}
		string /*NO &*/ code() const{return "\\ensuremath{"+code_+"}";}
		friend math operator<(math const& self, math const& other){
			return math(self.code_ + " < " + other.code_);
		}
		//friend math operator+(math const& sum1, math const& sum2){
		//	return math("{"+sum1.code_+"}+{"+sum2.code_+"}");
		//}
		friend math operator^(math const& base, math const& exponent){
			return math("{"+base.code_+"}^{"+exponent.code_+"}");
		}
		friend math sqrt(math const& base){
			return math("\\sqrt{"+base.code_+"}");
		}
		template<long N>
		friend math power(math const& base){
			return math("{"+base.code_+"}^{"+boost::lexical_cast<string>(N)+"}");
		}
		friend math log(math const& base){
			return math("\\log\\left(" + base.code_ +"\\right)");
		}
	};
	template<long N>
	math root(math const& base){
		return math("\\sqrt["+lexical_cast<string>(N)+"]{"+base.code_+"}");
	}
	math parenthesis(math const& contents){
		return math("{\\left("+contents.code_+"\\right)}");
	}
	struct plus_or_minus : math{
		plus_or_minus(string const& s) : math(s){}
		math operator^(math const& exponent) const{return math("\\left("+code_+"\\right)")^exponent;}
	};
	struct substraction : plus_or_minus{
		substraction(string s) : plus_or_minus(s){}
	};
	struct addition : plus_or_minus{
		addition(string s) : plus_or_minus(s){}
	};
	substraction operator-(math const& sub1, math const& sub2){
		return substraction(sub1.code_ + "-" + sub2.code_);
	}
	addition operator+(math const& sum1, math const& sum2){
		return addition(sum1.code_ + "+" + sum2.code_);
	}
	struct number : math{
		number(int const& i) : math(boost::lexical_cast<string>(i)){}
		number(double d) : math(lexical_cast<string>(d)){}
	};
	math operator-(math const& m){
		return math("-"+m.code_);
	}
	math exp(math const& exponent){
		return math("\\exp("+exponent.code_+")");
	}
	template<>
	math root<2>(math const& base){
		return math("\\sqrt{"+base.code_+"}");
	}
	struct product : math{
		product(math const& m) : math(m){}
		math operator^(math const& exponent) const{return math("\\left("+code_+"\\right)")^exponent;}
	};
	product operator*(math const& factor1, plus_or_minus const& factor2){
		return math("{"+factor1.code_+"}{"+parenthesis(factor2).code_+"}"); //assumme objects multiply like numbers, no symbol for product
	}
//	product operator*(math const& factor1, substraction const& factor2){
//		return math("{"+factor1.code_+"}{"+parenthesis(factor2).code_+"}"); //assumme objects multiply like numbers, no symbol for product
//	}
	product operator*(plus_or_minus const& factor1, math const& factor2){
		return product(parenthesis(factor1).code_+"{"+factor2.code_+"}");
	}
//	product operator*(addition const& factor1, math const& factor2){
//		return product(parenthesis(factor1).code_+"{"+factor2.code_+"}");
//	}
	product operator*(math const& factor1, math const& factor2){
		return math("{"+factor1.code_+"}{"+factor2.code_+"}"); //assumme objects multiply like numbers, no symbol for product
	}
	product operator*(number const& factor1, number const& factor2){
		return math("{"+factor1.code_+"}\\cdot{"+factor2.code_+"}");
	}
	struct sum : math{
		sum(math const& m) : math(m){}
		friend math operator*(sum const& s, math const& m){
			return math("("+s.code_+"){"+m.code_+"}");
		}
	};
	struct division : math{
		division(math const& m) : math(m){}
		math operator^(math const& exponent) const{return math("\\left("+code_+"\\right)")^exponent;}
	};
	division operator/(math const& factor1, math const& factor2){
		return math("\\frac{"+factor1.code_+"}{"+factor2.code_+"}"); //assumme objects multiply like numbers
	}

	//friend 
	ostream& operator<<(ostream& os, double const& d){
		(ostream::base&)(os) << math(d).code_; return os;
		//((base&)(os))<< d; return os;
	}
	template<class Unit> //friend 
	ostream& operator<<(ostream& os, boost::units::quantity<Unit> const& q){
		(ostream::base&)(os) << math(q).code(); return os;
	}
	template<class Dimension, class System>
	ostream& operator<<(ostream& os,boost::units::unit<Dimension, System> const& u){
		(ostream::base&)(os) << "\\ensuremath{" << u << "} "; return os;
	}

	struct environment {//: ostream{
	};
	template<class Environment>
	struct environment_category : environment{
		string label_;
		string options_;
		string arguments_;
		string postscript_;
		bool asterisk_;
		//virtual string preamble() const{return "";}
		environment_category() : asterisk_(false){}
		friend ostream& operator<<(ostream& os, Environment const& f){
			os.preamble()<< f.preamble();
			string environment_name = Environment::name();
			if(f.asterisk_) environment_name+="*";
			os << "\n";
			os << 
				"\\begin{" << environment_name << "}"<< ((f.options_=="")?(""):("["+ f.options_ +"]")) << ((f.arguments_=="")?(""):("{"+ f.arguments_ +"}")) << "\n" << 
					f.str() << "\n" << 
					((Environment::name()!="lstlisting")?("\\label{" + ((f.label_!="")?f.label_:boost::lexical_cast<string>(&f)) + "}\n"):("")) <<
				"\\end{" << environment_name << "}" << f.postscript_ << "\n";
			return os;
		}
		Environment& operator*(){asterisk_=true; return (Environment&)*this;}
		Environment& options(string const& opt){options_ = opt; return static_cast<Environment&>(*this);}
		Environment& operator[](string const& opt){return this->options(opt);}
	};
	template<class Figure>
	struct figure_category : environment_category<Figure>, ostream{
		//ostream{
		string caption;
		bool asterisk;
		figure_category() : asterisk(false){}
		using ostream::str;
		/*
		friend ostream& operator<<(ostream& os, figure_category<Figure> const& f){
			os.preamble()<< f.preamble();
			os << "\n\\begin{" << Figure::name() << (f.asterisk?"*":"") << "}\n"
				 << f.str() << "\n";
			if(not f.caption.empty()) os <<"\\caption{"<<f.caption<<"}\n";
			os << "\\end{" << Figure::name()   << (f.asterisk?"*":"") << "}\n";
			return os;
		}*/
	};
	struct figure : figure_category<figure>{
		figure& operator*(){asterisk=true; return *this;}
		static string name(){return "figure";}
	};
	template<class Command>
	struct command_category{
		string options_;
		friend ostream& operator<<(ostream& os, Command const& c){
			os.preamble()<< c.preamble();
			os << "\\" << Command::name() << ((c.options_!="")?("["+c.options_+"]"):("")) <<"{"
				 << c.str(); //<< "\n";
			os << "}\n";
			return os;
		}
		Command& options(string const& opt){options_=opt; return (Command&)*this;}
	};
	struct attachfile : ostream, command_category<attachfile>{
		attachfile(path const& p){
			assert(exists(p));
			(*this) << p.string();
		}
		//"\\attachfile[appearance=false, print=false, subject=Bibliography, description=BibTeX format, mimetype=text/plain]{" + bibfile + ".bib}";
		static string name(){return "attachfile";}
		protected:
		using ostream::str;
		friend ostream& operator<<(ostream& os, attachfile const& c);
	};
	struct embedfile : ostream, command_category<embedfile>{
		embedfile(path const& p){
			assert(exists(p));
			(*this) << p.string();
		}
		static string name(){return "embedfile";}
		protected:
		using ostream::str;
		friend ostream& operator<<(ostream& os, embedfile const& c);
	};
	struct marginnote : ostream, command_category<marginnote>{
		marginnote(){}
		marginnote(string contents){ (*this)<<contents;}
		static string name(){return "marginnote";}
		protected:
		using ostream::str;
		friend ostream& operator<<(ostream& os, marginnote const& c);
	};
	struct abstract: ostream, environment_category<abstract>{
		abstract(){}
		abstract(string contents){ (*this)<<contents;}
		static string name(){return "abstract";}
		using ostream::str;
	};
	struct equation : ostream, environment_category<equation>{
		//equation(string code){ (*this) << code;}
		equation(math const& m, std::string const& label=""){
			label_ = label;
			(*this) << m;
		}
		static string name(){return "equation";}
		using ostream::str;
	};
	struct tabular : ostream, environment_category<tabular>{
		static string name(){return "tabular";}
		template<class U1, class U2>
		tabular(std::vector<std::pair<boost::units::quantity<U1>, boost::units::quantity<U2> > > const& v){
			std::ostringstream coordinates_as_table;
			arguments_ = "D{.}{.}{-1}D{.}{.}{-1}"; //cc for center center, needs dcolumn "D{}{.}{5.4}D{.}{.}{5.4}"
			(*this) << "\\hline\\hline\n" << 
				// first row
				"\\multicolumn{1}{c}{" << boost::units::name<typename U1::dimension_type>() << "} & " << //the multicolumn spec is necessary when using dcolumn D{.}...
				"\\multicolumn{1}{c}{" << boost::units::name<typename U2::dimension_type>() << "} \\\\\n" << 
				// second row
				"\\multicolumn{1}{c}{" << "$(" << latex_string(U1()) << ")$} & " << 
				"\\multicolumn{1}{c}{" << "$(" << latex_string(U2()) << ")$}\\\\\n" << "\\hline\n";
			for(unsigned i = 0; i != v.size(); ++i){
				std::ostringstream oss_first; oss_first << std::setprecision(6); oss_first << v[i].first.value();
				std::ostringstream oss_second; oss_second << std::setprecision(6); oss_second << v[i].second.value();
				(*this) << "" << math(v[i].first.value()) << " & " << v[i].second.value() << "\\\\\n";
				coordinates_as_table<< v[i].first.value() << "\\t" << v[i].second.value() << "\\n";
			}
			(*this) << "\\hline\\hline";
			std::ostringstream postscript_oss;
			postscript_oss << "{\\pdfmarginpar[Note, Open=false, CA=0.05, Subject={" //<< "}, Caption={" //Author={au " << subject <<"} no implemented in pdfmarginpar 0.8 
						<< boost::units::name<typename U1::dimension_type>()<< " ["<< U1() <<"], "
						<< boost::units::name<typename U2::dimension_type>()<< " ["<< U2() <<"]  "
						<< "}]{"<< coordinates_as_table.str() << " }"
						<< "}\n";
			postscript_ = postscript_oss.str();
		}
		using ostream::str;
	};
	namespace tufte{
		struct fullwidth : 
			ostream, 
			environment_category<fullwidth>{
			fullwidth(){}
			fullwidth(string contents){ (*this) << contents;}
			static string name(){return "fullwidth";}
			friend struct environment_category<fullwidth>;
			using ostream::str;
		};
	}
	struct tikzpicture;
	struct tikz_object{
		virtual tikzpicture& operator>>(tikzpicture&) const=0;
	};
	struct tikzpicture : ostream{
		tikzpicture(){}
		tikzpicture(tikzpicture const& other) :  std::basic_ios<char, std::char_traits<char> >(), ostream(other.str()){}
		static string name(){return "tikzpicture";}
		//tikzpicture& operator<<(tikz_object const& t);//{assert(0); return *this;}
		friend tikzpicture& operator<<(tikzpicture& self, tikz_object const& t){
			t>>self;
			return self;
		}
		friend ostream& operator<<(ostream& os, tikzpicture const& tp){
			os.temp_files_.insert(tp.temp_files_.begin(), tp.temp_files_.end());
			os.preamble()<< tp.preamble();
			os << "\\begin{" << tikzpicture::name() << "}\n"
				 << tp.str() << "\n";
			os << "\\end{" << tikzpicture::name() << "}\n";
			return os;
		}
		template<class U1, class U2>
		tikzpicture& operator<<(std::vector<std::pair<boost::units::quantity<U1>, boost::units::quantity<U2> > > const&);
	};
	namespace tikz{
		typedef tikzpicture picture;
	}
	tikzpicture& operator<<(tikzpicture& tp, boost::multi_array<bool, 2> const& f){
		tp << "\\draw[help lines] ("<< -0.5 <<", "<< -0.5 <<") grid (" << f.shape()[0]-0.5 <<", " << f.shape()[1]-0.5 <<");";
		for(unsigned i = 0; i!=f.shape()[0]; ++i){
			for(unsigned j=0; j!=f.shape()[1]; ++j){
				if(f[i][j]==true){
					tp  << "\\draw[very thick, ->] (" << i+0. <<", "<<j-0.25 <<")--("<<i+0.<<", "<<j+0.25<<");";
				}else{
					tp  << "\\draw[very thick, <-] (" << i+0. <<", "<<j-0.25 <<")--("<<i+0.<<", "<<j+0.25<<");";
				}
			}
		}
		return tp;
	}
	struct marginfigure : figure_category<marginfigure>{
		static string name(){return "marginfigure";}
	};
	template<class CRTP>
	struct section_category{
		string title_;
		bool asterisk_;
		section_category(string title) : title_(title), asterisk_(false){}
		CRTP& operator*(){asterisk_=true; return static_cast<CRTP&>(*this);}
		friend ostream& operator<<(latex::ostream& os, CRTP const& self){
			os<<string("\n\n\\"+CRTP::name()+(self.asterisk_?string("*"):string(""))+"{");
			english title_english(self.title_);
			os<<title_english;// self.title_;
			return os<<string("}\n");
		}
	};
	struct section : public section_category<section>{
		explicit section(string title) : section_category<section>(title){}
		static string name(){return "section";}
	};
	struct subsection : section_category<subsection>{
		explicit subsection(string title) : section_category<subsection>(title){}
		static string name(){return "subsection";}
	};
	struct subsubsection : section_category<subsubsection>{ //not defined in tufte latex
		subsubsection(string title) : section_category<subsubsection>(title){}
		static string name(){return "subsubsection";}
	};
	struct newthought : section_category<newthought>{  // defined in tufte latex
		newthought(string title) : section_category<newthought>(title){}
		static string name(){return "newthought";} // needs a freaking double new-line before 
	};
	/// Automatic generation of (LaTeX) PGFplot's
	namespace listings{
		struct lstlisting : ostream, environment_category<lstlisting>{
			lstlisting(lstlisting const& other) : std::basic_ios<char, std::char_traits<char> >(), ostream(other.str()){}
			lstlisting(string contents){(*this) << contents;}
			lstlisting(char const contents[]){(*this) << contents;}
			static string name(){return "lstlisting";}
			using ostream::str;
		};
	}
	#define LATEX_LISTING_CODE(code) listings::lstlisting(#code); code
	namespace pgfplots{
		using boost::named;
		struct pair : std::pair<double, double>{
			pair(std::pair<double, double> const& p) : std::pair<double, double>(p){}
			pair(double const& d1, double const& d2) : std::pair<double, double>(d1, d2){}
			template<class Unit1, class Unit2> struct units;
		};
		template<class Unit1, class Unit2> 
		struct pair::units : pair{
			units(
				boost::units::quantity<Unit1> const& q1, 
				boost::units::quantity<Unit2> const& q2
			) : 
			pair(q1.value(), q2.value()){}
			template<class QA, class QB>
			units(
				//boost::units::quantity<UnitA> 
				QA const& qa, 
				//boost::units::quantity<UnitB> 
				QB const& qb
			)   :
				pair(
					boost::units::quantity<Unit1>(qa).value(), 
					boost::units::quantity<Unit2>(qb).value()
			){}
		};
		struct triple : boost::tuples::tuple<double, double, double>{
			triple(boost::tuples::tuple<double, double, double> const& t) : boost::tuples::tuple<double, double, double>(t){}
			triple(double const& d1, double const& d2, double const& d3) : boost::tuples::tuple<double, double, double>(d1, d2, d3){}
			template<class UX, class UY, class UZ> class units;
		};
		template<class UX, class UY, class UZ>
		struct triple::units : triple{
			units(
				boost::units::quantity<UX> const& q1, 
				boost::units::quantity<UY> const& q2, 
				boost::units::quantity<UZ> const& q3
			) : triple(
				q1.value(), 
				q2.value(), 
				q3.value()
			){} 
		};
		using std::vector; 
		typedef std::map<double, double> map;
		class axis_base;
		class axis;
		//class semilogxaxis; use "xmode = log, ymode = log" options to simplify
		struct plot : std::ostringstream{
			plot(){}
			plot(plot const& other) : 
				std::basic_ios<char, std::char_traits<char> >(), //solve warning in g++ 4.4.3 -Wextra
				std::ostringstream(other.str()){}
			virtual void operator>>(axis& a) const=0;
			virtual ~plot(){}
			template<class UnitX, class UnitY>
			class units;
		};
		struct axis_base : std::ostringstream, public tikz_object{
			//using std::ostringstream::operator<<;
			std::set<path> temp_files_;
			string options_;
			string coordinates_as_table_preamble_;
			axis_base(string opt="") : options_(opt){}
			axis_base& addlegendentry(string const& s){(*this) << "\\addlegendentry{" << s << "}\n"; return *this;}
			virtual tikzpicture& operator>>(tikzpicture& tp) const{
				//tp.temp_files_.insert(a.temp_files_.begin(), a.temp_files_.end());
				tp.usepackage("pgfplots");
				tp.preamble()<< "\\usepackage{pgfplots}\n\\usepgfplotslibrary{units}\n\\pgfplotsset{unit code/.code 2 args={{#1#2}}, unit marking pre={(}, unit marking post = {)} }\n"; // unit markings=parenthesis}\n";
				//os<<"\\tikzstyle{background rectangle}=[fill=white,opacity=0.1]\n"; //yellow
				tp<<
					"\\begin{axis}["<<this->options_<<"]\n"
						<<this->str()<<"\n"
					"\\end{axis}\n"
				;
				tp.preamble()<< "\\usetikzlibrary{backgrounds}\n";
				tp.preamble()<< "\\usetikzlibrary{patterns}\n";
				//tp.preamble()<< "\\usepgfplotslibrary{clickable}\n"; //needs insdljs.sty (AcroTeX)
				tp.preamble()<< "\\usepgfmodule{decorations}\n";
				return tp;
			}
			template<class UX, class UY>
			class units;
		};
		template<class UX, class UY>
		struct axis_base::units : public axis_base{
			units(string options="") : axis_base(options){}
		};
		template<class AxisT>
		struct axis_category : /*novirtual*/ axis_base{
			protected:
			axis_category(string options = "") : axis_base(options){}
			public:
			string code() const{
				return 
					"   \\begin{" + AxisT::name() + string("}[") +options_ + string("]\n")
					+		str() + string("\n"
					"   \\end{")  + AxisT::name() + string("}\n");
			}
			AxisT const& save(path p) const{
				clog << AxisT::name()<<"::save "<<p.string()<< endl;
				boost::filesystem::ofstream ofs(p);
				ofs << 
					"   \\begin{"+ AxisT::name() +"}["<<options_<<"]\n"
							<<str()<<"\n"
					"   \\end{"  + AxisT::name() +"}\n"
				;
				return static_cast<AxisT const&>(*this);
			}
			//friend ostream& operator<<(ostream& os, AxisT const& a){
			// bm1
			friend tikzpicture& operator<<(tikzpicture& os, AxisT const& a){ //axises can be put in a tikzpicture only
				os.temp_files_.insert(a.temp_files_.begin(), a.temp_files_.end());
				os.usepackage("pgfplots");
				os.preamble()<< "\\usepackage{pgfplots}\n\\usepgfplotslibrary{units}\n\\pgfplotsset{unit code/.code 2 args={{#1#2}}, unit marking pre={(}, unit marking post = {)} }\n"; // unit markings=parenthesis}\n";
				//os<<"\\tikzstyle{background rectangle}=[fill=white,opacity=0.1]\n"; //yellow
				os<<
					"\\begin{"+ AxisT::name()+"}["<<a.options_<<"]\n"
						<<a.str()<<"\n"
					"\\end{"+ AxisT::name() +"}\n"
				;
				os.preamble()<< "\\usetikzlibrary{backgrounds}\n";
				os.preamble()<< "\\usetikzlibrary{patterns}\n";
				//os.preamble()<< "\\usepgfplotslibrary{clickable}\n";
				os.preamble()<< "\\usepgfmodule{decorations}\n";
				return os;
			}
			template<class T>
			AxisT& operator<<(boost::named<T> const& c){
				static_cast<AxisT&>(*this) << "\\addlegendentry{"<<c.name()<<"}\n";
				c >> static_cast<AxisT&>(*this);
				return static_cast<AxisT&>(*this);
			}
			template<class UX, class UY>
			struct units;  // axis_category::units decl
		};
		class axis;
		class coordinates;
		class file; // pgfplots::file decl
		struct axis : 
			axis_category<axis>{
			typedef coordinates coordinates_type;
			struct options;  // axis::options decl
			static std::string name(){return "axis";}
			axis(string options = "") : axis_category<axis>(options){coordinates_as_table_preamble_ = /*"\\#"+*/ options +"\n";}// + options;}
			axis(options const& o);
			axis& operator()(std::map<string, coordinates> const& m);
			axis& operator<<(plot const& c){
				c >> *this; 
				return *this;
			}
			template<class T>
			axis& operator<<(boost::named<T> const& c){ //axis::operator<<
				(*this) << "\\addlegendentry{"<<c.name()<<"}\n";
				c >> *this;
				return *this;
			}
			template<class T>
			void operator<<(std::pair<string, T> const& p){
				std::string coordinates_as_table = p.second.str();
				boost::algorithm::replace_all(coordinates_as_table, "\n", "\\n");
				boost::algorithm::replace_all(coordinates_as_table, "(", "");
				boost::algorithm::replace_all(coordinates_as_table, ")", "");
				boost::algorithm::replace_all(coordinates_as_table, "%", "");
				boost::algorithm::replace_all(coordinates_as_table, ",", "\\t");
				std::string subject = p.first;
				(*this) << p.second;
				(*this) 
					<< "\\addlegendentry{"
					//<< "\\pdfmarkupcomment[author={author},subject={" << subject
					//<< "},color=red,opacity=0.01, markup=underline]{" << p.first << "}{" << coordinates_as_table << "} "
					<< p.first << " " 
					<< "\\pdfmarginpar[Note, Open=false, CA=0.05, Subject={" << subject <<" " //<< "}, Caption={"
					<< "}]{"<< coordinates_as_table << " }"
					<< "}";
			}
			template<class UnitX, class UnitY> struct units; // axis::units decl
		};
		struct axis::options : protected string, public std::map<string, string>{ // axis::options def
			std::map<double, string> ytick;
			std::map<double, string> xtick;
			options(){}
			options(string const& s) : string(s){}
			//options(char const* const& c) : string(c){}
			using std::map<string, string>::operator[];
			std::string code() const{
				string ret;
				ret += (string const&)(*this) + ", ";
				for(std::map<string, string>::const_iterator it = std::map<string, string>::begin(); it!= std::map<string, string>::end(); ++it){
					ret += it->first + " = " + it->second + ", ";
				}
				if(not xtick.empty()){
					string xtick_code; string xticklabels_code;
					for(std::map<double, string>::const_iterator it = xtick.begin(); it!=xtick.end(); ++it){
						if(it!=xtick.begin()){
							xtick_code+=", ";
							xticklabels_code+=", ";
						}
						xtick_code+=boost::lexical_cast<string>(it->first);
						xticklabels_code+="{" + it->second + "}";
					}
					ret+= "xtick = {"+xtick_code +"}, xticklabels = {"+ xticklabels_code +"}, ";
				}
				if(not ytick.empty()){
					string ytick_code; string yticklabels_code;
					for(std::map<double, string>::const_iterator it = ytick.begin(); it!=ytick.end(); ++it){
						if(it!=ytick.begin()){
							ytick_code+=", ";
							yticklabels_code+=", ";
						}
						ytick_code+=boost::lexical_cast<string>(it->first);
						yticklabels_code+="{" + it->second + "}";
					}
					ret+= "ytick = {"+ytick_code +"}, yticklabels = {"+ yticklabels_code +"}, ";
				}
				return ret;
			}
			options& operator|=(options const& other){
				string::operator+=(", ");
				string::operator+=(other);
				for(std::map<string, string>::const_iterator it = other.std::map<string, string>::begin(); it!= other.std::map<string, string>::end(); ++it){
					std::map<string, string>::insert(*it);
				}
				return *this;
			}
			friend options operator|(options const& op1, options const& op2){
				options ret(op1); return ret|=(op2);
			}
			friend options operator|(options const& op1, char const* const& op2){
				options ret(op1); return ret|=(options(op2));
			}
		};
		axis::axis(options const& o) : axis_category<axis>(o.code()){}
		static axis::options const no_markers("no markers");
		/// Semi logarithmic x axis  : use xmode = log
		template<class UX, class UY>
		struct plot::units : std::ostringstream{
			string options_;
			string postscript_;
			units(units const& other) :
				std::basic_ios<char, std::char_traits<char> >(), //solves warning g++ 4 -Wextra
				std::ostringstream(other.str()),
				options_(other.options_),
				postscript_(other.postscript_)
			{}
			units(string options="", string postscript="") : 
				options_(options), 
				postscript_(postscript)
			{
				//this->fill(' '); 
				this->precision(15); //don't use lower precision since it will present rugged texture in graphs with high mean mean
			}
			virtual void operator>>(axis_base::units<UX, UY>& a) const=0;
		};
		struct file : virtual plot{ // pgfplots::file def
			public:
			boost::filesystem::path p_;
			std::string options_;
			file(boost::filesystem::path const& p) :p_(p){}
			file(std::string options, boost::filesystem::path const& p) : p_(p), options_(options){}
			/*virtual*/ 
			void operator>>(axis& ax) const{
				ax << "\\addplot+[" << options_ <<"] file {" << p_.string() <<"} node[]{ \\attachfile[color=yellow,icon=Graph,print=false,zoom=false]{" << p_.string() << "} };";  
				return;
			}
		};
		struct coordinates : virtual plot{                                     // pgfplots::coordinates def
			typedef pair pair_type;
			string options_;
			string postscript_;
			std::ostringstream text_table_;
			unsigned size_;
			coordinates(coordinates const& other) : 
				std::basic_ios<char, std::char_traits<char> >(), //solves warning in g++ 4.4.3 -Wextra
				plot(other), 
				options_(other.options_), 
				postscript_(other.postscript_), 
				text_table_(other.text_table_.str()){}
			coordinates& operator=(coordinates const&){
				assert(0);
			}
			coordinates(string s="", string p="") : options_(s), postscript_(p){
				text_table_.fill (' '); 
				text_table_.precision(10);
			} //e.g. coordinates("", "\\closedcycle");
			coordinates& operator<<(pair const& p){
				text_table_ << p.first <<"	" <<p.second<<" \n";
				if((boost::math::isfinite)(p.first) and (boost::math::isfinite)(p.second)){
					*this<<"( "<<p.first<<", "<<p.second<<" )\n";
					++size_;
				}else{
					*this << "%( "<<p.first<<", "<<p.second<<" )  %discarted point due to not finite number condition \n";
				}
				return *this;
			}
			typedef std::map<double, double> map;
			coordinates& operator()(map const& m){
				for(map::const_iterator it=m.begin(); it!=m.end(); ++it){
					this->operator<<(pair(*it));
				}
				return *this;
			}
			virtual void operator>>(axis& a) const{
				std::string coordinates_as_table = this->str();
				boost::algorithm::replace_all(coordinates_as_table, "\n", "\\n");
				boost::algorithm::replace_all(coordinates_as_table, "(", "");
				boost::algorithm::replace_all(coordinates_as_table, ")", "");
				boost::algorithm::replace_all(coordinates_as_table, "%", "");
				boost::algorithm::replace_all(coordinates_as_table, ",", "\\t");
				a   << "%(not units) addplot\n"
					<< "\\addplot"
					<< (options_!=""?"+["+options_+"]":string())  
					<< " coordinates { \n" <<this->str()<<"} "
					<< postscript_ //<<" ; \n";
					<<" node{"
					//<<"\\pdfmarginpar[Name=Note, Open=false, CA=0.1, Subject={}, caption={}]{ " << coordinates_as_table << " }"
					<<"}"
					<<"; \n"
				;
			}
			template<class UnitX, class UnitY> struct units; // decl class axis::units
		};
		using boost::numeric::interval;
		using boost::optional;
		template<class UX, class UY> 
		struct axis::units :                                   //  axis::units def
			axis_base::units<UX,UY>{
			struct options : axis::options{                    // axis::options def
				options(string const& s="") : axis::options(s){}
				std::map<quantity<UX>, string> extra_x_tick;
				std::map<quantity<UY>, string> extra_y_tick;
				optional<interval<quantity<UX> > > restrict_x_to_domain;
				std::string code() const{
					//return (string)(*this);}
					//operator std::string() const{
					string ret;
					ret += (string const&)(*this);
					if(ret!="") ret+=", ";
					if(not extra_x_tick.empty()){
						string extra_x_ticks = "extra x ticks = {";
						string extra_x_tick_labels = "extra x tick labels = {";
						for(typename std::map<quantity<UX>, string>::const_iterator it=
								extra_x_tick.begin(); it!=
								extra_x_tick.end(); ++it
						){
							extra_x_ticks += lexical_cast<string>(it->first.value());
							extra_x_tick_labels += lexical_cast<string>(it->second);
							if(boost::next(it)!=extra_x_tick.end()){
								extra_x_ticks += ", ";
								extra_x_tick_labels +=", ";
							}
						}
						extra_x_ticks += "} ";
						extra_x_tick_labels +="} ";
						ret+=extra_x_ticks+", "+extra_x_tick_labels;
					}
					ret+=", ";
					if(not extra_y_tick.empty()){
						string extra_y_ticks = "extra y ticks = {";
						string extra_y_tick_labels = "extra y tick labels = {";
						for(typename std::map<quantity<UY>, string>::const_iterator it=
								extra_y_tick.begin(); it!=
								extra_y_tick.end(); ++it
						){
							extra_y_ticks += lexical_cast<string>(it->first.value());
							extra_y_tick_labels += lexical_cast<string>(it->second);
							if(boost::next(it)!=extra_y_tick.end()){
								extra_y_ticks += ", ";
								extra_y_tick_labels +=", ";
							}
						}
						extra_y_ticks += "} ";
						extra_y_tick_labels +="} ";
						ret+=extra_y_ticks+", "+extra_y_tick_labels;
					}
					if(restrict_x_to_domain) ret+=", restrict x to domain = "+lexical_cast<string>(lower(*restrict_x_to_domain).value())+":"+lexical_cast<string>(upper(*restrict_x_to_domain).value());
					return ret;
				}
			};
			typedef coordinates::units<UX, UY> coordinates_type; //cannot name it "coordinates"
			typedef pair::units<UX, UY> pair_type;
			static std::string name(){return "axis";}
			units(string opts) //this can benefit from delegating constructors
				:
				axis_base::units<UX, UY>(
					string(latex_string(UX())!=""? //handle dimensionless units
						"x unit = {"
							"\\href{http://www.wolframalpha.com/input/?i=(" + name_string(UX()) + ")}{"
								+ latex_string(UX()) + "}}, ":"")
								+ string(latex_string(UY())!=""?
						"y unit = {" 
							//"\\href{http://www.wolframalpha.com/input/?i=(" + name_string(UY()) + ")}"
							"{" + latex_string(UY()) + "}}, ":"")+
						"xlabel = {" + boost::units::name<typename UX::dimension_type>() + "}, "
						"ylabel = {" + boost::units::name<typename UY::dimension_type>() + "}  " 
							+ ((options(opts).code().empty()?"":", ") + options(opts).code())
			){}
			units(options const& opt = options()) 
				: 
				axis_base::units<UX, UY>(
					string(latex_string(UX())!=""? //handle dimensionless units
						"x unit = {"
							"\\href{http://www.wolframalpha.com/input/?i=(" + name_string(UX()) + ")}{"
								+ latex_string(UX()) + "}}, ":"")
								+ string(latex_string(UY())!=""?
					"y unit = {" 
					"\\href{http://www.wolframalpha.com/input/?i=(" + name_string(UY()) + ")}{"
					+ latex_string(UY()) + "}}, ":"")+
				"xlabel = {" + boost::units::name<typename UX::dimension_type>() + "}, "
				"ylabel = {" + boost::units::name<typename UY::dimension_type>() + "}  " 
				+ ((opt.code().empty()?"":", ") + opt.code())
			){}
			units& operator<<(plot::units<UX, UY> const& c){
				c >> static_cast<units&>(*this); return static_cast<units&>(*this);
			}
			template<class T>
			units& operator<<(boost::named<T> const& c){ //axis::units::operator<<
				(*this) << "\\addlegendentry{"<<c.name()<<"}\n";
				c >> *this;
				return *this;
			}
			void save(path p) const{
				boost::filesystem::ofstream ofs(p);
				ofs<< "%\\usepackage{pgfplots}\n";
				ofs << 
					"\\begin{tikzpicture}\n"
					"   \\begin{axis}["<<options_<<"]\n"
							<< this->str()<<"\n"
					"   \\end{axis}\n"
					"\\end{tikzpicture}\n";
				return; // static_cast<AxisT const&>(*this);
			}
			template<class T>
			void operator<<(std::pair<char const*, T> const& p){
				(*this)<< std::make_pair(std::string(p.first), p.second);
			}
			template<class T>
			void operator<<(std::pair<string, T> const& p){
				std::string coordinates_as_table = p.second.str();
				boost::algorithm::replace_all(coordinates_as_table, "\n", "\\n");
				boost::algorithm::replace_all(coordinates_as_table, "(", "");
				boost::algorithm::replace_all(coordinates_as_table, ")", "");
				boost::algorithm::replace_all(coordinates_as_table, "%", "");
				boost::algorithm::replace_all(coordinates_as_table, ",", "\\t");
				std::string subject = p.first;
				//author={author}
				const_cast<T&>(p.second).node_at_end_ += "\\pdfmarkupcomment[author={dataset},color=Green, opacity=0.05,markup=Squiggly]{ }{"+subject + "\\n" + boost::lexical_cast<std::string>(p.second.size_) + " data points" +"}\n";
				//if(p.second.str()=="") (*this)<<"\\begin{comment} %due to lack of points (empty string check)\n";
				if(p.second.str()!=""){
					(*this) << p.second;
					(*this) 
						<< "\\addlegendentry{"
						<<"\\hypertarget{" << &p.second << "}{ }"
						//<< "\\pdfmarkupcomment[author={Donald Duck},subject={" << subject << " "
						//<< boost::units::name<typename UX::dimension_type>()<< " ["<< UX() <<"], "
						//<< boost::units::name<typename UY::dimension_type>()<< " ["<< UY() <<"]  "
						//<< "},color=red,opacity=0.1, markup=Underline]{" << p.first << "}{" << coordinates_as_table << "} "
						<< p.first << " " 
						<< "\\pdfmarginpar[Note, Open=false, CA=0.05, Subject={" << subject << " " //<< "}, Caption={" //Author={au " << subject <<"} no implemented in pdfmarginpar 0.8 
						<< boost::units::name<typename UX::dimension_type>()<< " ["<< UX() <<"], "
						<< boost::units::name<typename UY::dimension_type>()<< " ["<< UY() <<"]  "
						<< "}]{"<< coordinates_as_table << " }"
						<< "}\n"
					;
				}else{
					(*this) << "% addplot and addlegendentry omitted for dataset named " << p.first << " due to lack of points (empty string check)\n"; //if not commented it creates a missmatch
					std::clog << "% addplot and addlegendentry omitted for dataset named " << p.first << " due to lack of points (empty string check)\n";
				}
				//if(p.second.str()=="") (*this)<<"\\end{comment} %due to lack of points (empty string check)\n";
			}
			/*
			friend ostream& operator<<(ostream& os, units const& a){
				os.temp_files_.insert(a.temp_files_.begin(), a.temp_files_.end());
				os.usepackage("pgfplots");
				os.preamble()<< "\\usepackage{pgfplots}\n\\usepgfplotslibrary{units}\n\\pgfplotsset{unit code/.code 2 args={{#1#2}}, unit marking pre={(}, unit marking post = {)} }\n"; //unit markings=parenthesis}\n";
				os<<"\\tikzstyle{background rectangle}=[fill=white,opacity=0.1]\n"; //yellow
				os<<"	\\begin{"+ units::name()+"}["<<a.options_<<"]\n"
				  <<a.str()<<"\n";
				os<<"	\\end{"+ units::name() +"}      \n";
				os.preamble()<< "\\usetikzlibrary{backgrounds}\n";
				os.preamble()<< "\\usetikzlibrary{patterns}\n";
				os.preamble()<< "\\usepgfplotslibrary{clickable}\n";
				os.preamble()<< "\\usepgfmodule{decorations}\n";
				return os;
			}*/
		};
		using boost::units::quantity;
		template<class UX, class UY> 
		struct coordinates::units :                         // coordinates::units def
			plot::units<UX, UY>{
			string node_at_end_;
			unsigned size_;
			typedef pair::units<UX, UY> pair_type;
			std::ostringstream text_table_;
			units(units const& other) : 
				std::basic_ios<char, std::char_traits<char> >(), //solves g++ 4 -Wextra
				plot::units<UX, UY>(other),
				text_table_(other.text_table_.str())
			{}
			units(string options="", string postscript="") : 
				plot::units<UX, UY>(options, postscript)
			{ //options_(options), postscript_(postscript){ //coordinates(s, p){
				text_table_ << "# [" << UX() <<"], [" << UY() << "]\n";
			}
			template<  //or be more specific template<class system,  > ... (unit<system>...)
				class U1, //=UX, 
				class U2 //=UY  //function template default allowed in c++0x only
			>
			units& operator()(
				std::vector<std::pair<quantity<U1>, quantity<U2> > > const& v
			){
				for(unsigned idx=0; idx!=v.size(); ++idx){
					(*this) << pair_type(
						v[idx].first, 
						v[idx].second
					);
				}
				return *this;
			}
			template<class U1, class U2>
			units& operator()(std::map<quantity<U1>, quantity<U2> > const& m){
				for(typename std::map<quantity<U1>, quantity<U2> >::const_iterator it=m.begin(); it!=m.end(); ++it){
					(*this) << pair_type(it->first, it->second);
				}
				return *this;
			}
			virtual void operator>>(axis_base&) const{
				assert(0);
			}
			units& operator<<(pair::units<UX, UY> const& p){
				text_table_ << p.first << "\t" << p.second<<" \n";
				if((boost::math::isfinite)(p.first) and (boost::math::isfinite)(p.second)){
					*this<<"( "<<p.first<<", "<<p.second<<" )\n";
					++size_;
				}else{
					*this << "%( "<<p.first<<", "<<p.second<<" )\n"; 
				}
				return *this;
			}
			virtual void operator>>(axis_base::units<UX, UY>& a) const{
				std::string coordinates_as_table = this->str();
				boost::algorithm::replace_all(coordinates_as_table, "\n", "\\n");
				boost::algorithm::replace_all(coordinates_as_table, "(", "");
				boost::algorithm::replace_all(coordinates_as_table, ")", "");
				boost::algorithm::replace_all(coordinates_as_table, "%", "");
				boost::algorithm::replace_all(coordinates_as_table, ",", "\\t");
				a
					<< "\\addplot"
					<<(this->options_!=""?"+["+this->options_+"]":string(""))  
					<<" coordinates { \n"<<this->str()<<"} "
					<<this->postscript_ //<<" ; \n";
					<<" node{"
					<< node_at_end_ <<" " //<<"\\hyperlink{" << this << "}{ }" //space
					//<<"\\pdfmarginpar[Name=Note, Open=false, CA=0.01, Subject={}]{ " << coordinates_as_table << " }" //, caption={} caption works with pdfmarginpar 0.9
					<<"}"
					<<"; \n"
				;
			}
		};
		/* doesnt work
		template<class UX, class UY> 
		axis_base::units<UX, UY>& operator<<(
			axis_base::units<UX, UY>& axu, 
			std::pair<
				std::string, 
				pgfplots::axis::units<UX, UY>
			> const& p
		){
			return axu;
		}*/
		axis& axis::operator()(std::map<string, coordinates> const& m){
			(*this)<<"%axis::operator()(std::map<string, coordinates> const& m)\n";
			for(std::map<string, coordinates>::const_iterator it=m.begin(); it!=m.end(); ++it){
				std::string coordinates_as_table = it->second.str();
				boost::algorithm::replace_all(coordinates_as_table, "\n", "\\n");
				boost::algorithm::replace_all(coordinates_as_table, "(", "");
				boost::algorithm::replace_all(coordinates_as_table, ")", "");
				boost::algorithm::replace_all(coordinates_as_table, "%", "");
				boost::algorithm::replace_all(coordinates_as_table, ",", "\\t");
				//coordinates_as_table = this->coordinates_as_table_preamble_ + coordinates_as_table;
				//std::clog << "preamble is " << boost::algorithm::replace_all_copy(coordinates_as_table_preamble_, "\\rm", "") << std::endl;
				(*this) << "%no units addplot \n"
					<< "\\addplot"
					<<(it->second.options_!=""?"+["+it->second.options_+"]":string())  
					<<" coordinates { \n"<<it->second.str()<<"} "
					<<it->second.postscript_
					<<" node{"
					<<"}"
					<<"; \n"
				;
				(*this)
					<<"\\addlegendentry{"<<it->first
					//<<"\\pdfmarginpar[Note, Open=false, CA=0.1, Subject={" << it->first << "}, caption={"<< 
					//	boost::algorithm::replace_all_copy(boost::algorithm::replace_all_copy(coordinates_as_table_preamble_, ",", "\\n"), "\\rm ","")
					//<<"}]{ " << coordinates_as_table << " }"
					<<"}\n"
				;
			}
			return *this;
		}
		class coordinates3 : public plot {// {, public std::ostringstream{
			string preamble_;
			public:
			coordinates3(string p="") : preamble_(p){}
			coordinates3& operator<<(
				triple const& t
			){
				*this<<"( "
					<<boost::tuples::get<0>(t)<<", "
					<<boost::tuples::get<1>(t)<<", "
					<<boost::tuples::get<2>(t)
					<<") \n"; 
				return *this;
			}
			void operator>>(axis/*_base*/& a) const{
				a << "\\addplot3+["<<preamble_<<"] coordinates {\n" <<this->str()<<"};\n";
			}
			template<class UnitX, class UnitY, class UnitZ> struct units;
		};
		template<class UX, class UY, class UZ> 
		struct coordinates3::units : 
			public plot::units<UX, UY> {
			units(string options="") : plot::units<UX, UY>(options){}
			units& operator<<(triple::units<UX, UY, UZ> const& t){
				*this <<"( "
					<<boost::tuples::get<0>(t)<<", "
					<<boost::tuples::get<1>(t)<<", "
					<<boost::tuples::get<2>(t)
					<<" )\n";
				return *this;
			}
			virtual void operator>>(axis_base::units<UX, UY>& a) const{
				a << "\\addplot3+["<<this->options_<<"] coordinates {\n" << this->str()<<"};\n";
			}
		};
	}//ns pgfplots
	namespace pgf{
		namespace plots = pgfplots;
	}

	template<class U1, class U2>
	tikzpicture& tikzpicture::operator<<(std::vector<std::pair<boost::units::quantity<U1>, boost::units::quantity<U2> > > const& vp){
		typedef pgfplots::axis::units<U1, U2> axis; axis ax;
		typedef typename axis::coordinates_type coordinates; coordinates coords;
		ax << coords(vp);
		(*this) << ax;
		return *this;
	}
	template<class U1, class U2, class UA, class UB>
	pgfplots::axis::units<U1, U2>&  operator<<(pgfplots::axis::units<U1, U2>& ax, std::vector<std::pair<boost::units::quantity<UA>, boost::units::quantity<UB> > > const& vp){
		typename pgfplots::axis::units<U1, U2>::coordinates_type coords;
		coords(vp);
		ax << coords;
		return ax;
	}

}//ns latex

/*
template<class Functor, typename IntervalUnit>
latex::tikz::picture plot(
	Functor const& f, 
	boost::units::quantity<IntervalUnit> const& x1, 
	boost::units::quantity<IntervalUnit> const& x2
){
	using namespace latex;
	typedef boost::units::quantity<IntervalUnit> domain_type;
	typedef decltype(f(x1)) result_type; 
	typedef typename pgf::plots::axis::units<IntervalUnit, typename result_type::unit_type > axis;
	axis ax;
	typedef typename axis::coordinates_type coordinates; 
	coordinates coords("no markers")
	;
	for(
		boost::units::quantity<IntervalUnit> x=x1; 
		x<x2; 
		x+=(x2-x1)/100.
	){
		coords << typename coordinates::pair_type(x, f(x));
	}
	ax << coords;
	tikz::picture tkzp; 
	tkzp << ax;
	return tkzp;
}
*/
#endif
#ifdef _TEST_LATEX_HPP
#define ECHO(s) "\\hspace{-0.5cm}" << (#s) << std::endl << (s)
#include<iostream>
//#include<boost/assign.hpp>
//#include "boost/units/systems/atomic.hpp"
//#include "boost/units/systems/atomic/si_conversion.hpp"

#include <boost/mpl/string.hpp>

using std::cout; using std::endl;
using namespace boost::units;
using namespace latex;
latex::ostream lo("latex.pdf");
int main(){
	lo.usepackage("lipsum");
	lo.author = "Alfredo Correa";
	lo.title = "The C++ LaTeX Library";
	lo << maketitle;
	abstract abs("Abstract \\lipsum[1]");
	lo << abs;
	lo << section("Getting started")
		<< "To start writting a LaTeX file that ultimatelly will be compiled as a PDF, you need to construct an stream object with the name of the destination PDF. "
		<< "Title information of the document can be added optionally. "
		<< "\\begin{lstlisting}[frame=tb,language=C++]\n"
		<< "#include \"latex.hpp\"\n"
		<< " ... \n"
		<< "latex::ostream lo(\"latex.pdf\"); \n"
		<< "lo.author = \"A. U. Thor\"; \n "
		<< "lo.title = \"The C\\!++ \\LaTeX{} Library\";\n"
		<< "lo << latex::maketitle; \n"
		<< "\\end{lstlisting}\n"
	;
	lo << "The LaTeX output stream can accept any explicit command "
		<< "\\begin{lstlisting}[frame=tb,language=C++]\n"
		<< "lo << \"\\section{Getting Started}\\n\";"
		<< "\\end{lstlisting}\n"
		<< "or alternatively any of the predefined latex commands and environements\n"
		<< "\\begin{lstlisting}[frame=tb,language=C++]\n"
		<< "lo << latex::section(\"Getting Started\");"
		<< "\\end{lstlisting}\n"
		<< "The ubiquitous latex namespace specifier can be eliminated by declaring the use of the namespace "
		<< "\\begin{lstlisting}[frame=tb,language=C++]\n"
		<< "using namespace latex;\n" 
		<< "lo << section(\"Getting Started\");"
		<< "\\end{lstlisting}\n"
		<< "It is called output stream (ostream) because it works similarly to an standard output stream, "
			"that is, the PDF is only produced (updated) if the flush and endl sent to the stream or if "
			"the stream is destroyed (or goes out of scope). For example"
		<< "\\begin{lstlisting}\n"
		<< "lo << \"Last line\" << endl;\n"
		<< "\\end{lstlisting}\n"
		<< "will add such line to the document, and end of line character and then will flush the result and generate the PDF. "
		<< "If endl or flush is used too frequently the resulting program will be slow because it will compile the PDF each time. ";
		lo << "To end a line/paragraph, the string \"\\\\n\\\\n\", newline, or par can be used instead. " << par << 
		"Long paragraphs can be define by using the default C++ literal string concatenation. "
		"The latex control object par can be used to break paragraphs." << newline << "And newline can be used to break lines."
	;
	lo << section("Output of C++ results") <<
		"Things get interesting we start using the output for dynamical C++ results in a program."
		"For example, the following C++ code" <<
		listings::lstlisting(
			"for(double x = 0.; x<5.; ++x) lo << \"$\" << x << \"^2 = \" << x*x << \"$, \";"
		)
		<< "Outputs the following line" << par;
		for(double x = 0; x < 5; ++x) lo << "$" << x << "^2 = " << x*x << "$, ";
	listings::lstlisting l = (
		"latex::equation quadratic(\"x^2\");\n"
		"lo << quadratic << \"as seen in previous equation \"<< latex::ref(&quadratic) <<\". \";\n"
	);
	latex::equation quadratic("x^2");
	lo << section("Cross References") <<
		"Crossed reference can be managed by explicitly giving \\verb+\\label+/\\verb+\\ref+ pairs as usual in \\LaTeX{}, "
		"alternatively pointer reference to the objects can be used as references. "
		"This code " << l << "generates this text" << 
		quadratic << "as seen in previous equation " << latex::ref(&quadratic) <<". "
	;
	lo << LATEX_LISTING_CODE(
		for(double x = 0.; x<5.; ++x){
			lo <<  x << ", ";
		}
	);
	lo << par;
	lo << quantity<si::length>::from_value(23) << newline;
	lo << quantity<si::dimensionless>::from_value(12) << newline;
	tikz::picture p;
	boost::multi_array<bool, 2> ma(boost::extents[4][4]);
	bool cma[4][4] = {{true,false,true,false},{true,false,true,true},{false,false,true,false},{true,true,false,true}};
	ma.assign(cma, cma+16); 
	for(unsigned i=0; i!=4; ++i){
		for(unsigned j=0; j!=4; ++j){
		ma[i][j] = cma[i][j];
		cout << ma[i][j] << endl;
		}
	}
	p << ma;
	lo << ECHO( number(4)*number(5) = number(20) ); 
	lo << p;
	{
		pgfplots::axis a("small, title = a");
		boost::named<pgfplots::coordinates> c("c", pgfplots::coordinates()); 
		c 
			<< pgfplots::pair(1.,2.)
			<< pgfplots::pair(3.,4.);
		a << c;
		tikzpicture tkzp; tkzp << a;
		lo << tkzp;
	}
	{
		pgfplots::axis a("small, title = {Plotting from file}");
		pgfplots::file f("no markers, blue", "/home/correaa/prj/alf/xyz.gr");
		a << f;
		tikz::picture p; p << a;
		lo << p;
	}
	{
		lo << section("Spell Checker")
			<< "The library supports spell checking which is implemented by marking mispelled words in the resulting pdf. In the following Wikipedia quote some words have been mispelled. (In order to see the corrections and suggestions you need to open the PDF with Acrobat Reader; to avoid printing the corrections from Acrobat Reader in the pring box select to print Document only, instead of Document and Markup)" << par;
		lo << 
			english("Shakespeare was vorn and raized in \\emph{Stratford-upon-Avon}. At the age of 18, he narried Anne Hathaway\\cite{caca}, with whom he had three children: Susanna, and twins Hamnet and Judith. Between 1585 and 1592, he began a succesful career in London as an acto, writer, and part owner of a playing company called the Lord Chamberlain's Men, later known as the King's Men. He appears to have retirred to Stratford around 1613, where he died three years later. Few records of Shakespeare's private life survive, and there has been considerable speculacion about such maters as his physical appearance, sexuality, religious beliefs, and whether the works attributed to him were writen by others.");
	}
	/*
	{
		using namespace pgfplots;
		//typedef pgfplots::axis::units<si::time, si::length> axis;
		axis::units<si::time, si::length> a("title = b, legend pos = south east");
		//boost::named<axis::coordinates_type> c(axis::coordinates_type());
		boost::named<coordinates::units<si::time, si::length> > c("a", coordinates::units<si::time, si::length>()); //(
			//std::string("a"), 
			//axis::units<si::time, si::length>::coordinates_type());
		c << pair::units<si::time, si::length>(1.*si::second, 2.*si::meter);
		c << pair::units<si::time, si::length>(2.*si::second, 3.*si::meter);
		//axis::pair_type(1.*si::second, 2.*si::meter);
		//c << axis::pair_type(2.*si::second, 3.*si::meter);
		(a << c);
		tikz::picture p; p << a;
		lo << p;
	}*/
	//lo << ECHO("hola") << std::endl << ECHO(1+3) << std::endl;
	#if 0
	lo <<"chau"<<std::endl;
	lo << latex::par << 0.4*atomic::hartree/atomic::hbar << latex::par;
	std::map<double, double> m=boost::assign::map_list_of
		(2,-2.8559703)
		(3,-3.5301677)
		(4,-4.3050655)
		(5,-5.1413136)
		(6,-6.0322865)
		(7,-6.9675052)
		(8,-7.9377747);
	std::map<double, double> m2=boost::assign::map_list_of
		(2,-3.8559703)
		(3,-2.5301677)
		(4,-6.3050655)
		(5,-2.1413136)
		(6,-1.0322865)
		(7,-8.9675052)
		(8,-1.9377747);
	latex::pgfplots::axis a;
	lo << "\\pgfplotsset{"
		"legend style={"
		"font=\\footnotesize,"
		"fill=white, opacity = 0.5"
		"}}\n";
	a << latex::pgfplots::named_coordinates( latex::pgfplots::coordinates()<< m, "m");
	//a << boost::named<latex::pgfplots::coordinates&>(latex::pgfplots::coordinates()<< m, "m");
	a << (latex::pgfplots::coordinates()<< m2);
	//	a 
	//	<< (latex::pgfplots::coordinates() << m) 
	//	<< (latex::pgfplots::coordinates() << m2); 
	lo << a << "\n\n hola \n\n" << a;
	latex::pgfplots::semilogxaxis b;
	b << (latex::pgfplots::coordinates() << m);
	lo << b;
	{
		using namespace boost::units;
		using namespace atomic;
		cout << latex_string(bohr*bohr*bohr*bohr*bohr*bohr) <<endl;
	}
	{
		quantity<divides<nonsi::electron_volt_unit, nonsi::angstrom_unit>::type> f(3.*si::joules/si::meter);
		quantity<si::force> fsi(f);
		latex::pgfplots::pair::units<si::length, si::force> p(
			1.*si::meter,
			//3.*si::joules/si::meter
			f
		);
		quantity<si::length> qs(1.*atomic::bohr);
	}
	{
		typedef make_scaled_unit<si::time, scale<10, static_rational<-15> > >::type fs_unit;
		quantity<fs_unit> t;
		quantity<divide_typeof_helper<si::length, fs_unit>::type > v;
		cout << t << endl; 
		cout << v << endl; 
	}
	#endif
	cout << boost::lexical_cast<std::string, float>(0.05) << endl;
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

