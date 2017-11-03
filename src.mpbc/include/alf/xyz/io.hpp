#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -Wall `#-Wfatal-errors` -Wextra -Wno-unused-variable .$0.cpp -I$HOME/usr/include -I$HOME/prj -L$HOME/usr/lib -lboost_filesystem -lboost_system -D_TEST_XYZ_IO_HPP -o .$0.cpp.x && ./.$0.cpp.x $1 $2 $3 $4 $5 $6 $7 $8 $9
exit
#endif
#ifndef XYZ_IO_HPP
#include "../xyz.hpp"
#include "../asymptote.hpp"
namespace xyz{
	asymptote::ostream& operator<<(asymptote::ostream& ao, frame const& f){
		using namespace asymptote;
		using namespace graph3;
		double scale = 100.;
		ao.import("three");
		unsigned n = f.size();
		if(not f.empty()){
			triple t1 = {{(*f)[0][0]*scale, (*f)[0][1]*scale, (*f)[0][2]*scale}};
			ao << "begingroup3(\"crystal\");\n";
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[0][0]*scale, (*f)[0][1]*scale, (*f)[0][2]*scale}}),"gray+linewidth(5), currentlight");
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[1][0]*scale, (*f)[1][1]*scale, (*f)[1][2]*scale}}),"gray+linewidth(5), currentlight");
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[2][0]*scale, (*f)[2][1]*scale, (*f)[2][2]*scale}}),"gray+linewidth(5), currentlight");
			ao << "endgroup3();\n";
		}
		#ifdef XYZ_IO_NO_ATOMS
		{
			ao << "begingroup3(\"atoms\");\n";
				for(unsigned i=0; i!=n; ++i){
					using namespace asymptote::three;
					surface sphere("scale3(70)*unitsphere"); //pm
					triple position={{ f[i][0]/nonsi::angstrom*scale, f[i][1]/nonsi::angstrom*scale, f[i][2]/nonsi::angstrom*scale }};
					ao.draw(shift(position)*sphere, "yellow, name = \""+ f[i] +"\""); //render(compression=Low, merge=true)
				}
			ao << "endgroup3();\n";
		}
		#endif
		ao << "begingroup3(\"bonds\");\n";
		for(unsigned i=0; i!=n; ++i){
			triple position={{ f[i][0]/nonsi::angstrom*scale, f[i][1]/nonsi::angstrom*scale, f[i][2]/nonsi::angstrom*scale }};
			for(unsigned j=i+1; j!=n; ++j){
				triple position_j={{ f[j][0]/nonsi::angstrom*scale, f[j][1]/nonsi::angstrom*scale, f[j][2]/nonsi::angstrom*scale }};
				array<double, 3> d = {{position_j[0] - position[0], position_j[1] - position[1], position_j[2] - position[2]}};
				double s = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
				if(sqrt(s)<2.*scale){
					ao.draw(connect(position,position_j),"yellow+linewidth(40), currentlight"); //, render(compression=High, merge=true)
				}
			}
		}
		ao << "endgroup3();\n";
		return ao;
	}
	void export_asymptote(frame const& f, boost::filesystem::path p){
		using namespace asymptote;
		using namespace graph3;
		ostream ao(p);
		ao.import("three");
		unsigned n = f.size();
		if(not f.empty()){
			triple t1 = {{(*f)[0][0]*100., (*f)[0][1]*100., (*f)[0][2]*100.}};
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[0][0]*100., (*f)[0][1]*100., (*f)[0][2]*100.}}),"gray+linewidth(5), currentlight");
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[1][0]*100., (*f)[1][1]*100., (*f)[1][2]*100.}}),"gray+linewidth(5), currentlight");
			ao.draw(connect((triple){{0,0,0}},(triple){{(*f)[2][0]*100., (*f)[2][1]*100., (*f)[2][2]*100.}}),"gray+linewidth(5), currentlight");
		}
		for(unsigned i=0; i!=n; ++i){
			using namespace asymptote::three;
			surface sphere("scale3(70)*unitsphere"); //pm
			triple position={{ f[i][0]/nonsi::angstrom*100., f[i][1]/nonsi::angstrom*100., f[i][2]/nonsi::angstrom*100. }};
			ao.draw(shift(position)*sphere, "yellow");
			for(unsigned j=i+1; j!=n; ++j){
				triple position_j={{ f[j][0]/nonsi::angstrom*100., f[j][1]/nonsi::angstrom*100., f[j][2]/nonsi::angstrom*100. }};
				array<double, 3> d = {{position_j[0] - position[0], position_j[1] - position[1], position_j[2] - position[2]}};
				double s = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
				if(sqrt(s)<200.){
					ao.draw(connect(position,position_j),"yellow+linewidth(40), currentlight");
				}
			}
		}
	}
}
#endif
#ifdef _TEST_XYZ_IO_HPP
int main(){
	xyz::ifstream xyzi("beta-boron.xyz");
	xyz::frame f; xyzi >> f;
	//asymptote::ostream ao("beta-boron.pdf");
	export_asymptote(f, "beta-boron.pdf");
	return 0;
}
#endif

