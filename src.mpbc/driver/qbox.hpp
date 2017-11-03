#ifndef QBOX_HPP_
#define QBOX_HPP_
#define CLOG(expr) (std::clog<<(#expr)<<" = "<<(expr)<<" ("<<__FILE__<<":"<<__LINE__<<": "<<__FUNCTION__<<")"<<std::endl)
#define BOOST_ENABLE_ASSERT_HANDLER
#include "alf/boost/multi_array_exception.hpp"

//#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
//#include<boost/lexical_cast.hpp>
#include<boost/string_cast.hpp>
#include<boost/numeric/conversion/cast.hpp>
#include<boost/math/constants/constants.hpp>
#ifndef _NO_BOOST_MPI
#include<boost/mpi.hpp>
#endif
#include<boost/multi_array.hpp>
#include<boost/noncopyable.hpp>
#include<boost/optional.hpp>

//#include<boost/units/static_rational.hpp>
#include<boost/units/quantity.hpp>
#include"boost/units/systems/atomic.hpp"
#include"boost/units/systems/atomic/si_conversion.hpp"
#include<boost/units/systems/si/io.hpp> //beautify units output
#include<boost/units/systems/si/pressure.hpp> //si::pascal
#include<boost/units/systems/si/time.hpp> //si::second
#include<boost/units/make_scaled_unit.hpp>
#include<boost/units/systems/si/codata/physico-chemical_constants.hpp> //m_u (dalton unit, http://en.wikipedia.org/wiki/Atomic_mass_unit) and k_B (boltzman constant)

//#include<boost/units/base_units/metric/angstrom.hpp> //for angstrom is xyz export

#ifndef _NO_BOOST_MPI
#include"boost/mpi/ostream.hpp"
#endif
#include<boost/random.hpp>

#include "alf/boost/array_io.hpp"
#include "alf/fftw3.hpp"
#include "alf/hdf5.hpp"
#include "alf/mod.hpp"
#include "alf/xyz.hpp"

#include "Basis.h"
#include "FourierTransform.h"
#include "BOSampleStepper.h"
#include "Sample.h"
#include "SampleReader.h"
#include "SampleWriter.h"
#include "SlaterDet.h"
#include "Species.h"
#include "SpeciesReader.h"
#include "UnitCell.h"
#include "Wavefunction.h"
#include "IonicStepper.h"

#include "boost/units/systems/atomic.hpp"
#include "boost/units/systems/atomic/si_conversion.hpp"
//units
namespace qbox{
using namespace boost::units;
typedef quantity<atomic::energy > energy;
typedef quantity<atomic::time   > time;
typedef quantity<si::temperature> temperature;
typedef quantity<atomic::volume > volume;
typedef boost::units::make_scaled_unit<boost::units::si::pressure, boost::units::scale<10, boost::units::static_rational<9> > >::type gpa_unit;
typedef boost::units::quantity<gpa_unit> pressure;
static const gpa_unit gigapascal, GPa; //maybe define conversion factor to pascal or atomic::pressure
using nonatomic::rydberg;
using nonatomic::Ry;
using atomic::hartree;
using atomic::bohr;
}

#include "alf/linear_space.hpp"
#include "field.hpp"

#include<iostream>
#include<iomanip>
#include<stdexcept>



#include "qbox/vector.hpp"
#include "qbox/cell.hpp"

#include "UserInterface.h"

//this is to process commands from files
//??namespace{
using std::cout;
using std::string;
using std::min;
using std::max;
using std::ios;
using std::setprecision;
using std::ostringstream;
using std::setw;
//using std::vector;//can't
#include "AtomsDyn.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "ChargeMixCoeff.h"
#include "ChargeMixRcut.h"
#include "Debug.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecuts.h"
#include "Emass.h"
#include "ExtStress.h"
#include "FermiTemp.h"
#include "Dt.h"
#include "Nempty.h"
#include "NetCharge.h"
#include "Nrowmax.h"
#include "RefCell.h"
#include "Stress.h"
#include "Thermostat.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "Xc.h"

#include "AngleCmd.h"
#include "AtomCmd.h"
#include "ComputeMLWFCmd.h"
#include "ConstraintCmd.h"
#include "DistanceCmd.h"
#include "FoldInWsCmd.h"
#include "HelpCmd.h"
#include "KpointCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "MoveCmd.h"
#include "PrintCmd.h"
#include "QuitCmd.h"
#include "RandomizeWfCmd.h"
#include "ResetVcmCmd.h"
#include "RunCmd.h"
#include "SaveCmd.h"
#include "SetCmd.h"
#include "SpeciesCmd.h"
#include "StatusCmd.h"
#include "StrainCmd.h"
#include "TorsionCmd.h"
#include "AtomsDyn.h"


//?? }//ns unnamed

namespace qbox {
//using std::cout;
using std::endl;
using boost::numeric_cast;
using boost::string_cast;
//using boost::units::static_rational;
using boost::multi_array_ref;
typedef std::complex<double> complex;
//static energy const rydberg(0.5 * atomic::hartree);
class energies {
public:
	energies(EnergyFunctional const& ef) :
		ekin_ (ef.ekin()    *hartree), //provides units to qbox numbers
		econf_(ef.econf()  *hartree),
		ehart_(ef.ehart()  *hartree),
		ecoul_(ef.ecoul()  *hartree),
		exc_(ef.esr()      *hartree),
		eself_(ef.eself()  *hartree),
		ets_(ef.ets()      *hartree),
		etotal_(ef.etotal()*hartree) {}
	energy
			ekin_,
	       econf_,
	       eps_,
	       enl_,
	       ehart_,
	       ecoul_,
	       exc_,
	       esr_,
	       eself_,
	       ets_,
	       etotal_;
	energy ekin() const {
		return ekin_;
	}
	energy econf() const {
		return econf_;
	}
	energy eps() const {
		return eps_;
	}
	energy enl() const {
		return enl_;
	}
	energy ehart() const {
		return ehart_;
	}
	energy ecoul() const {
		return ecoul_;
	}
	energy exc() const {
		return exc_;
	}
	energy esr() const {
		return esr_;
	}
	energy eself() const {
		return eself_;
	}
	energy ets() const {
		return ets_;
	}
	energy etotal() const {
		return etotal_;
	}
};
class state : public energies {
public:
	boost::optional<pressure> p_;
	state(BOSampleStepper const& ss) :
		energies(ss.ef()),
		p_(ss.s_.ctrl.stress == "ON",
		   quantity<gpa_unit>((
							   (ss.sigma_eks[0] + ss.sigma_eks[1] + ss.sigma_eks[2]) +
							   (ss.sigma_kin[0] + ss.sigma_kin[1] + ss.sigma_kin[2])
							   ) / 3.0 * atomic::E_h / pow<3>(atomic::a_0)
							  )){}
	pressure p() const {
		return *p_;   //will raise an exception if pressure was not calculated!
	}
};
class session : boost::noncopyable {
private:
	Context ctxt;
public:
	Sample* s;
private:
	std::ostream& out;
	mutable std::ostringstream history_;
	UserInterface ui;
	void init() {
		//init variables as in original code
		s->ctrl.atoms_dyn = "LOCKED";
		//s->ctrl.cell    = 0;
		s->ctrl.cell_dyn  = "LOCKED";
		s->ctrl.cell_lock = "OFF";
		s->ctrl.charge_mix_coeff = 0.5;
		s->ctrl.charge_mix_rcut = 10.0;
		s->ctrl.cell_mass  = 10000.0;
		s->ctrl.debug = "OFF";
		//s->ctrl.ecut = 0.0;
		s->ctrl.ecutprec   = 0.0;
		s->ctrl.ecuts      = 0.0;
		s->ctrl.emass      = 0.0;
		for(int i = 0; i < 6; i++) s->ctrl.ext_stress[i] = 0.0;
		s->ctrl.fermi_temp = 0.0;
		s->ctrl.dt         = 3.0;
		//s->ctrl.nempty = 0;
		//s->ctrl.net_charge = 0.0;
		s->wf.set_nrowmax(1);
		s->wf.update_occ(0.0);
		//s->ctrl.ref_cell = 0;
		s->ctrl.stress     = "OFF";
		s->ctrl.thermostat = "OFF";
		s->ctrl.th_temp    = 0.0;
		s->ctrl.th_width   = 100.;
		s->ctrl.wf_diag    = "F";
		s->ctrl.wf_dyn     = "SD";
		s->ctrl.xc         = "LDA";

		//init user interface for interpreting input files
		ui.addCmd(new AngleCmd(s));
		ui.addCmd(new AtomCmd(s));
		ui.addCmd(new ComputeMLWFCmd(s));
		ui.addCmd(new ConstraintCmd(s));
		ui.addCmd(new DistanceCmd(s));
		ui.addCmd(new FoldInWsCmd(s));
		ui.addCmd(new HelpCmd(s));
		ui.addCmd(new KpointCmd(s));
		ui.addCmd(new ListAtomsCmd(s));
		ui.addCmd(new ListSpeciesCmd(s));
		ui.addCmd(new LoadCmd(s));
		ui.addCmd(new MoveCmd(s));
		ui.addCmd(new PrintCmd(s));
		ui.addCmd(new QuitCmd(s));
		ui.addCmd(new RandomizeWfCmd(s));
		ui.addCmd(new ResetVcmCmd(s));
		ui.addCmd(new RunCmd(s));
		ui.addCmd(new SaveCmd(s));
		ui.addCmd(new SetCmd(s));
		ui.addCmd(new SpeciesCmd(s));
		ui.addCmd(new StatusCmd(s));
		ui.addCmd(new StrainCmd(s));
		ui.addCmd(new TorsionCmd(s));
		ui.addVar(new AtomsDyn(s));
		ui.addVar(new Cell(s));
		ui.addVar(new CellDyn(s));
		ui.addVar(new CellLock(s));
		ui.addVar(new ChargeMixCoeff(s));
		ui.addVar(new ChargeMixRcut(s));
		ui.addVar(new CellMass(s));
		ui.addVar(new Debug(s));
		ui.addVar(new Ecut(s));
		ui.addVar(new Ecutprec(s));
		ui.addVar(new Ecuts(s));
		ui.addVar(new Emass(s));
		ui.addVar(new ExtStress(s));
		ui.addVar(new FermiTemp(s));
		ui.addVar(new Dt(s));
		ui.addVar(new Nempty(s));
		ui.addVar(new NetCharge(s));
		ui.addVar(new Nrowmax(s));
		ui.addVar(new RefCell(s));
		ui.addVar(new Stress(s));
		ui.addVar(new Thermostat(s));
		ui.addVar(new ThTemp(s));
		ui.addVar(new ThTime(s));
		ui.addVar(new ThWidth(s));
		ui.addVar(new WfDiag(s));
		ui.addVar(new WfDyn(s));
		ui.addVar(new Xc(s));
	}
public:
	session(std::ostream& os = std::cout) : // ?move to front
		ctxt(),
		s(new Sample(ctxt)), //replace by auto_ptr?
		out(os) {
		init();
	}
	~session() {
		delete s;
	}
	unsigned natoms() const {
		unsigned count = 0;
		for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
			count += s->atoms.atom_list[is].size();
		}
		return count;
	}
	AtomSet const& atoms() const {
		return s->atoms;
	}
	AtomSet& atoms() {
		return s->atoms;
	}
	Atom const& atom(unsigned idx) const {
		assert(idx < natoms());
		unsigned count = 0;
		for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
			for(unsigned ia = 0; ia != s->atoms.atom_list[is].size(); ++ia) {
				if(count == idx) {
					return *(s->atoms.atom_list[is][ia]);
				}
				count++;
			}
		}
		return *(s->atoms.atom_list[0][0]);
	}
	session& atom(
	    std::string name, 
		std::string species,
	    vector const& position,
	    vector velocity = vector(0, 0, 0)
	) {
		history_ 
			<< "atom " << name << " " << species << " " 
			<< position << " ";
		if(velocity!=vector(0,0,0))
			history_ << velocity;
		history_<<std::endl;
		Atom* a = new Atom(name, species, position, velocity);
		const int atoms_nel_before = s->atoms.nel();
		if(!(s->atoms.addAtom(a))) {
			delete a;
			throw std::logic_error("could not add atom `" + name + "'");
		}
		const int atoms_nel_after = s->atoms.nel();
		const int delta_nel = atoms_nel_after - atoms_nel_before;
		const int wf_nel = s->wf.nel();
		s->wf.set_nel(wf_nel + delta_nel);
		s->wf.update_occ(0.0);
		if(s->wfv != 0) {
			s->wfv->set_nel(wf_nel + delta_nel);
			s->wfv->clear();
		}
		return *this;
	}
	session& import_xyz(boost::filesystem::path const& p, std::map<std::string, boost::filesystem::path> const& species);
	template<class OStream>
	session const& export_xyz(OStream& oss, std::string comment) const { //add mpbc option
		//typedef atomic::angstrom_base_unit::unit_type angstrom_unit; //this should be in xyz_format::angstrom
		//static const nonsi::angstrom_unit angstrom; //BOOST_UNITS_STATIC_CONSTANT(angstrom, angstrom_unit);
		using nonsi::angstrom_unit;
		using nonsi::angstrom;
		{
			unsigned count = 0;
			for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
				count += s->atoms.atom_list[is].size();
			}
			oss << count << '\n' << comment 
				<< " # vol = " << get_cell().volume()*pow<3>(atomic::a_0)
				<< " # coordinates in angstroms" << '\n';
		}
		unsigned line = 0;
		for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
			for(unsigned ia = 0; ia != s->atoms.atom_list[is].size(); ++ia) {
				oss << "  "
				    << s->atoms.atom_list[is][ia]->name()[0] << " "
				    << quantity<angstrom_unit>(s->atoms.atom_list[is][ia]->position().x * bohr).value() << " "
				    << quantity<angstrom_unit>(s->atoms.atom_list[is][ia]->position().y * bohr).value() << " "
				    << quantity<angstrom_unit>(s->atoms.atom_list[is][ia]->position().z * bohr).value() << " ";
				// for compatibility with mathematica disable the following
				if(line == 0) {
					oss
					        << " crystal_vector 1 " << s->atoms.cell().a(0)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom)
					        << " crystal_vector 2 " << s->atoms.cell().a(1)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom)
					        << " crystal_vector 3 " << s->atoms.cell().a(2)*(double)(quantity<angstrom_unit>(1.*atomic::bohr) / angstrom);
				}
				if(line == 1) {
					oss << " bbox_xyz " //assumes orthogonal cell
					    << " 0. " << s->atoms.cell().a(0).x* (double)(quantity<angstrom_unit>(1.*bohr) / angstrom)
					    << " 0. " << s->atoms.cell().a(1).y* (double)(quantity<angstrom_unit>(1.*bohr) / angstrom)
					    << " 0. " << s->atoms.cell().a(2).z* (double)(quantity<angstrom_unit>(1.*bohr) / angstrom);
				}
				oss << '\n';
				++line;
			}
		}
		oss << std::flush;
		return *this;
	}
	session const& export_xyz(boost::filesystem::path const& p, std::string comment="# -1") const;
	void kpoint_add(vector const& v, double weight) {
		history_ << "kpoint add " << v << " " << weight << endl;
		s->wf.add_kpoint(v, weight);
	}
	//using namespace boost::filesystem;
	//atoms positions as qbox commands
	session const& export_sys(boost::filesystem::path const& p, std::string comment="# comment, atomic units positions/velocity") const{
		history_ << "# export sys " << p << endl;
		boost::filesystem::ofstream ofs(p);
		ofs << "set cell ";
		for(unsigned i=0; i!=3; ++i){
			for(unsigned j=0; j!=3; ++j){
				ofs<< s->atoms.cell().a(i)[j]<<" ";
			}
		}
		ofs << "\n";
		//ofs << "set ref_cell ";
		//for(unsigned i=0; i!=3; ++i){
		//	for(unsigned j=0; j!=3; ++j){
		//		ofs<< s->ctrl.ref_cell.a(i)[j]<<" ";
		//	}
		//}
		for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
			ofs << "species " << s->atoms.species_list[is]->name() << " " << s->atoms.species_list[is]->uri() << '\n';
		}
		for(unsigned is = 0; is != s->atoms.species_list.size(); ++is) {
			for(unsigned ia = 0; ia != s->atoms.atom_list[is].size(); ++ia) {
				ofs << "atom " 
				    << s->atoms.atom_list[is][ia]->name() << " "
				    << s->atoms.atom_list[is][ia]->species() <<" "
				    << s->atoms.atom_list[is][ia]->position().x<< " "
				    << s->atoms.atom_list[is][ia]->position().y<< " "
				    << s->atoms.atom_list[is][ia]->position().z<< "   "
				    << s->atoms.atom_list[is][ia]->velocity().x<< " "
				    << s->atoms.atom_list[is][ia]->velocity().y<< " "
				    << s->atoms.atom_list[is][ia]->velocity().z<< "\n"
				;
			}
		}
		return *this;
	}
	void fold_in_ws(){
		s->atoms.fold_in_ws();
	}
	void kpoint_delete(vector const& v) {
		history_ << "kpoint delete " << v << endl;
		s->wf.del_kpoint(v);
	}
	void randomize_wf(double amp = 0.02) {
		history_ << "randomize_wf" << endl;
		s->wf.randomize(amp);
	}
	session& set_cell(cell const& c) {
		history_ << "set cell " << c.a(0) << "  " << c.a(1) << "  " << c.a(2) << endl;
		s->wf.resize(c, s->wf.refcell(), s->wf.ecut());
		if(s->wfv != 0) {
			s->wfv->resize(c, s->wf.refcell(), s->wf.ecut());
			s->wfv->clear();
		}
		s->atoms.set_cell(c.a(0), c.a(1), c.a(2));
		return *this;
	}
	cell get_cell() const {
		return s->wf.cell();   //can't have the name cell()
	}
	session& set_ecut(energy ecut) {
		return set_ecut( quantity<atomic::dimensionless>(ecut / rydberg) );
	}
	session& set_ecut(double ecut /*rydberg*/) { // ? make private
		history_ << "set ecut " << ecut << endl;
		s->wf.resize(0.5 * ecut); //internally uses atomic units (hartree) internally
		if(s->wfv) {
			s->wfv->resize(0.5 * ecut);
			s->wfv->clear();
		}
		return *this;
	}
	session& set_ecutprec(energy ecutprec) {
		return set_ecutprec(quantity<atomic::dimensionless>(ecutprec / rydberg));
	}
	session& set_ecutprec(double v /*rydberg*/) {
		history_ << "set ecutprec " << v << endl;
		s->ctrl.ecutprec = 0.5 * v; //internally uses atomic units (hartree) internally
		return *this;
	}
	session& set_ecuts(energy ecuts) {
		return set_ecuts( quantity<atomic::dimensionless>(ecuts / rydberg));
	}
	session& set_ecuts(double v /*rydberg*/) {
		history_ << "set ecuts " << v << endl;
		assert(v > 0);
		s->ctrl.ecuts = 0.5 * v;
		//assert(s->ctrl.ecuts<ecut_);
		return *this;
	}
	session& set_nrowmax(unsigned v) {
		history_ << "set nrowmax " << v << endl;
		s->wf.set_nrowmax(v);
		s->wf.update_occ(0.0);
		if(s->wfv != 0) {
			s->wfv->set_nrowmax(v);
			s->wfv->clear();
		}
		return *this;
	}
	session& set_xc(std::string xc = "LDA") {
		history_ << "set xc " << xc << endl;
		s->ctrl.xc = xc;
		return *this;
	}
	session& set_wf_dyn(std::string dyn) { //give default?
		history_ << "set wf_dyn " << dyn << endl;
		s->ctrl.wf_dyn = dyn;
		return *this;
	}
	session& set_atoms_dyn(std::string dyn = "LOCKED") {
		assert(dyn == "LOCKED" or dyn == "SD" or dyn == "SDA" or dyn == "CG" or dyn == "MD");
		history_ << "set atoms_dyn " << dyn << endl;
		s->ctrl.atoms_dyn = dyn;
		return *this;
	}
	session& set_cell_dyn(std::string dyn = "LOCKED" ) {
		history_ << "set cell_dyn " << dyn << endl;
		s->ctrl.cell_dyn = dyn;
		return *this;
	}
	session& set_charge_mix_coeff(double const& d = 0.5) {
		assert(d > 0.0 and d <= 1.0);
		history_ << "set carge_mix_coeff " << d << endl;
		s->ctrl.charge_mix_coeff = d;
		return *this;
	}
	double get_charge_mix_coeff() const {
		return s->ctrl.charge_mix_coeff;
	}
	void set_charge_mix_rcut(double const& d) {
		history_ << "set carge_mix_rcut " << d << endl;
		s->ctrl.charge_mix_rcut = d;
	}
	session& set_dt(double d = 3.0) { // ? make private
		history_ << "set dt " << d << endl;
		s->ctrl.dt = d;
		return *this;
	}
	session& set_dt(time d){
		return set_dt(d/atomic::time_unit);
	}
	qbox::time dt() const {
		return qbox::time(s->ctrl.dt * atomic::time_unit);
	}
	session& set_fermi_temp(temperature const& t) {
		return set_fermi_temp(t / si::kelvin);
	}
	session& set_fermi_temp(double d = 0.0) {
		s->ctrl.fermi_temp = d;
		return *this;
	}
	session& set_nempty(unsigned v = 0) {
		history_ << "set nempty " << v << endl;
		s->wf.set_nempty(v);
		if(s->wfv != 0) {
			s->wfv->set_nempty(v);
		}
		return *this;
	}
	unsigned nempty() const {
		return s->wf.nempty();
	}
	session& set_stress(bool v) {
		history_ << "set stress ";
		if(v) {
			history_ << "ON";
			s->ctrl.stress = "ON";
		} else {
			history_ << "OFF";
			s->ctrl.stress = "OFF";
		}
		history_ << std::endl;
		return *this;
	}
	session& reset_vcm() {
		history_ << "reset_vcm" << std::endl;
		s->atoms.reset_vcm();
		return *this;
	}
	session& set_thermostat(string v = "OFF") {
		history_ << "set thermostat " << v << std::endl;
		assert((v == "SCALING") or(v == "ANDERSEN") or(v == "LOWE") or(v == "OFF"));
		s->ctrl.thermostat = v;
		return *this;
	}
	session& set_th_temp(double t) {
		assert(t >= 0);
		history_ << "set th_temp " << t << " # kelvin" << std::endl;
		s->ctrl.th_temp = t;
		return *this;
	}
	session& set_th_temp(temperature T) {
		return set_th_temp(T / si::kelvin);
	}
	session& set_th_time(double tau = 5000. /*au of time = 120 fs*/) {
		assert(tau >= 0);
		history_ << "set th_time " << tau << " # atomic unit of time" << std::endl;
		s->ctrl.th_time = tau;
		return *this;
	}
	session& set_th_time(time tau) {
		return set_th_time(tau / atomic::time_unit);
	}
	session& set_th_width(double const T = 100. /*kelvin*/) {
		assert(T >= 0);
		history_ << "set th_width " << T << " # kelvin" << std::endl;
		s->ctrl.th_time = T;
		return *this;
	}
	session& set_th_width(temperature T) {
		return set_th_width(T / si::kelvin);
	}
	session& species(
	    std::string name,
	    boost::filesystem::path const& p
	) {
		history_ << "species " << name << " " << p << endl;
		using boost::filesystem::exists;
		using std::runtime_error;
		if(not exists("species.xsd")) {
			throw runtime_error("file `species.xsd' does not exists");
		}
		if(not exists(p)) {
			throw runtime_error("file `" + p.string() + "' does not exists");
		}
		SpeciesReader sp_reader(s->ctxt_);
		Species* sp = new Species(s->ctxt_, name);
		sp_reader.readSpecies(*sp, p.string());
		sp_reader.bcastSpecies(*sp);
		s->atoms.addSpecies(sp, name);
		return *this;
	}
	session& status() {
		s->wf.info(cout, "wf");
		if(s->wfv != 0) {
			s->wfv->info(cout, "wfv");
		}
		cout << "<vcm> " << s->atoms.vcm() << " </vcm>" << endl;
		return *this;
	}
	session& load(boost::filesystem::path const& p) try {
		history_ << "load " << p << endl;
		using boost::filesystem::exists;
		using std::runtime_error;
		if(not exists(p)) {
			throw runtime_error("sample file '" + p.string() + "' does not exists to load");
		}
		if(not exists("sample.xsd")) {
			throw runtime_error("schema file 'sample.xsd' not found");
		}
		if(not exists("species.xsd")) {
			throw runtime_error("schema file 'species.xsd' not found");
		}
		SampleReader reader(s->ctxt_);
		reader.readSample(*s, p.string(), false);
		s->ctxt_.barrier();
		return *this;
	} catch(std::exception& e){
		throw std::runtime_error("can not read '" + p.string() + "' because "+e.what());
	} catch(... /*std::exception& e*/) {
		throw std::runtime_error("can not read '" + p.string() + "' due to unknown reason" /*, "+e.what()*/);
	}
	session& process_commands(boost::filesystem::path const& p) {
		history_ << p.string() << " #external script";
		if(not boost::filesystem::exists(p)) {
			history_ << ", warning, not found!";
		}
		history_ << endl;
		boost::filesystem::ifstream ifs(p);
		ui.processCmds(ifs, const_cast<char*>(("[" + p.string() + "]").c_str()), true);
		return *this;
	}
	session& save(boost::filesystem::path const& p, std::string const& description="saved with qbox driver") {
		history_ << "save " << p << endl;
		{
			bool base64 = true;
			bool atomsonly = false;
			bool serial = false;
			bool save_wfv = true;
			SampleWriter(s->ctxt_).writeSample(
											   *s, p.string(),
											   description, base64, atomsonly, serial, save_wfv
											   );
		}
		return *this;
	}
	session& save_atomsonly(boost::filesystem::path const& p, std::string const& description="saved with qbox driver only atoms"){
		clog << "save_atomsonly saves the file but it can't be loaded with the load command" << endl;
		history_<<"save -atomsonly "<<p<<endl;
		{
			bool base64 = true;
			bool atomsonly = true; //other setting are irrelevant if atomsonly=true;
			bool serial = true; //false;
			bool save_wfv = false; //true;
			SampleWriter(s->ctxt_).writeSample(
											   *s, p.string(),
											   description, base64, atomsonly, serial, save_wfv
											   );
		}
		return *this;
	}
	session& export_history(boost::filesystem::path const& p) {
		boost::filesystem::ofstream ofs(p);
		std::ofstream(p.string().c_str());
		ofs << history_.str();
		return *this;
	}
	qbox::state	run(int niter, int nitscf = 1, int nite = 1) { //qbox::run object?
		//std::clog << "begin run" << std::endl;
		assert(s->atoms.species_list.size()>0 and s->atoms.atom_list[0].size()>0); //check that we are not running with no atoms (typicall cause of corrupted restart file)
		history_ << "run " << niter << " " << nitscf << " " << nite << endl;
		BOSampleStepper bo(*s, nitscf, nite, out);
		s->wf.info(out, "wavefunction");
		//std::clog << "run0: ptr " << bo.ionic_stepper << std::endl;
		bo.step(niter);
		//std::clog << "run1: ptr " << bo.ionic_stepper << std::endl;
		return qbox::state(bo);
	}
	Wavefunction const& wavefunction() const {
		return s->wf;
	}
	double external_magnetic_field() const {
		int flux =
#ifndef _MPBC
		    0
#else
		    1
#endif
		    ;
		cell c = get_cell();
		double area = length(c.a(1) ^ c.a(2));
		std::clog << "the area is " << area << " (Bohr(a.u.)^2)" << std::endl;
		//double const pi = 3.141592;
		using boost::math::constants::pi;
		double b = 2 * pi<double>() * flux / area;
		std::clog << "magnetic field in au " << b << " = " << b * 235051.7 << " Tesla" << std::endl;
		return b;
	}
};

std::vector<double> const& occupations(Wavefunction const& Wf_, unsigned ispin = 0, unsigned ikp = 0) {
	return Wf_.sd(ispin, ikp)->occ_;
}
unsigned nempty(Wavefunction const& Wf_) {
	return Wf_.nempty();
}

using boost::const_multi_array_ref;
using boost::multi_array;
using boost::extents;
using boost::array;
using fftw3::plan;
using fftw3::backward;

class euclidean {
	qbox::vector v_;
	euclidean(qbox::vector const& v) : v_(v) {}
	operator qbox::vector() {
		return v_;
	}
};

class grid { //fft grid
public:
	grid(cell const& uc, boost::array<unsigned, 3> const& s) : domain_(uc), shape_(s) {}
	UnitCell domain_;
	boost::array<unsigned, 3> shape_;
	friend std::ostream& operator<<(std::ostream& os, grid const& self) {
		return os << "grid(" << self.domain_ << ", (array<unsigned, 3>){{" << self.shape_[0] << "," << self.shape_[1] << "," << self.shape_[2] << "}})";
	}
};

using boost::multi_array;
using boost::const_multi_array_ref;
using boost::extents;
using boost::array;
template<typename T, typename D = multi_array<T, 3u> >  //needs multi_array.hpp
class field {
public:
	field(grid const& g) :
		domain_(g.domain_),
		data_(boost::extents[g.shape_[0]][g.shape_[1]][g.shape_[2]]) {
	}
	T operator()(boost::array<double, 3> const& a) const {
		return operator()(vector(a[0], a[1], a[2]));
	}
	T operator()(vector const& v) const {
		//trilinear interpolation
		boost::array<int, 3> idxi;
		boost::array<double, 3> idxd;
		for(size_t i = 0; i != 3; ++i) {
			using boost::math::constants::pi;
			boost::tie(idxi[i], idxd[i]) = alf::modf((domain_.b(i) * v / (2 * pi<double>())) * data_.shape()[i]);
		}
		int const xf = idxi[0] % (data_.shape()[0]); //qbox::cell a and b are normalized to 2pi
		int const yf = idxi[1] % (data_.shape()[0]);
		int const zf = idxi[2] % (data_.shape()[0]);
		int const  xc = ((xf + 1) % (data_.shape()[0]));
		int const  yc = ((yf + 1) % (data_.shape()[1]));
		int const  zc = ((zf + 1) % (data_.shape()[2]));
		double const& xd = idxd[0];
		double const& yd = idxd[1];
		double const& zd = idxd[2];
		//from http://en.wikipedia.org/wiki/Trilinear_interpolation#Method
		T i1 = data_[xf][yf][zf] * (1 - zd) + data_[xf][yf][zc] * zd;
		T i2 = data_[xf][yc][zf] * (1 - zd) + data_[xf][yc][zc] * zd;
		T j1 = data_[xc][yf][zf] * (1 - zd) + data_[xc][yf][zc] * zd;
		T j2 = data_[xc][yc][zf] * (1 - zd) + data_[xc][yc][zc] * zd;
		T w1 = i1 * (1 - yd) + i2 * yd;
		T w2 = j1 * (1 - yd) + j2 * yd;
		return (w1 * (1 - xd) + w2 * xd);
	}
	T max() const { //T must be comparable
		double ret = data_.origin()[0];
		for(unsigned i = 0; i != data_.num_elements(); ++i) {
			if(ret < data_.origin()[i]) {
				ret = data_.origin()[i];
			}
		}
		return ret;
	}
	UnitCell const& domain() const {
		return domain_;
	}
	UnitCell domain_;
	D data_;
};
template<typename T>
T max(field<T> const& f) {
	return f.max();
};
//general case
template<class T>
hdf5::file const& operator<<(
    hdf5::file const& f,
    hdf5::nvp<field<T> > const& nf) {
	boost::const_multi_array_ref<double, 2> cell_as_ma(nf.const_value().domain().amat(), boost::extents[3][3]);
	boost::const_multi_array_ref<T, 3> den(nf.const_value().data_);
	f["cell"] << cell_as_ma;
	f["data"] << den;
	return f;
}
//specialized for vectors fields
template<class T>
hdf5::file const& operator<<(hdf5::file const& f, hdf5::nvp<field<boost::array<T, 3> > > const& nf) {
	const_multi_array_ref<double, 2> cell_as_ma(nf.const_value().domain().amat(), boost::extents[3][3]);
	const_multi_array_ref<double, 4> data((double*)nf.const_value().data_.origin(), boost::extents[nf.const_value().data_.shape()[0]][nf.const_value().data_.shape()[1]][nf.const_value().data_.shape()[2]][3]);
	f["cell"] << cell_as_ma;
	f["data"] << data;
	return f;
}
template<>
hdf5::file const& operator<<(hdf5::file const& f, hdf5::nvp<field<vector> > const& nf) {
	boost::const_multi_array_ref<double, 2> cell_as_ma(nf.const_value().domain().amat(), boost::extents[3][3]);
	boost::const_multi_array_ref<double, 4> data((double*)nf.const_value().data_.origin(), boost::extents[nf.const_value().data_.shape()[0]][nf.const_value().data_.shape()[1]][nf.const_value().data_.shape()[2]][3]);
	f["cell"] << cell_as_ma;
	f["data"] << data;
	return f;
}
hdf5::file const& operator<<(hdf5::file const& f, hdf5::nvp<field<complex> > const& nf) {
	const_multi_array_ref<double, 2> cell_as_ma(nf.const_value().domain().amat(), boost::extents[3][3]);
	const_multi_array_ref<double, 4> data((double*)nf.const_value().data_.origin(), boost::extents[nf.const_value().data_.shape()[0]][nf.const_value().data_.shape()[1]][nf.const_value().data_.shape()[2]][2]);
	f["cell"] << cell_as_ma;
	f["data"] << data;
	return f;
}
grid get_grid(Wavefunction const& wf_) {
	return grid(
	    wf_.cell(),
		(array<unsigned, 3>) {{
			wf_.sd(0, 0)->basis().np(0),
			wf_.sd(0, 0)->basis().np(1),
			wf_.sd(0, 0)->basis().np(2)
		}}
	);
}
template<
	class RationalFlux, 
	class MultiArray = const_multi_array_ref<std::complex<double>, 3> 
>
class wf_ref {
protected:
	grid grid_;
	MultiArray wfXYZ_ref_;
	wf_ref(
		grid const& g, 
		MultiArray const& data
	) : grid_(g), 
		wfXYZ_ref_(data) {}
public:
	wf_ref(Wavefunction const& Wf_, unsigned istate, unsigned ispin = 0, unsigned ikp = 0) :
		grid_(get_grid(Wf_)),
		wfXYZ_ref_(Wf_.sd(ispin, ikp)->c().cvalptr(istate* Wf_.sd(ispin, ikp)->c().mloc()),
		           extents
		           [Wf_.sd(ispin, ikp)->basis().np(0)]
		           [Wf_.sd(ispin, ikp)->basis().np(1)]
		           [Wf_.sd(ispin, ikp)->basis().np(2)]
		          ) {
	}
	cell domain() const { return grid_.domain_; }
	field<complex> wavefunction() const;
	field<double> density() const;
	field<vector> current() const;
	//field<array<double,3> > current() const;
};

template<class RationalFlux>
class wf :
	protected multi_array<complex, 3>,
	public wf_ref<RationalFlux> {
public:
	wf(grid const& g) :
		multi_array<complex, 3>(extents[g.shape_[0]][g.shape_[1]][g.shape_[2]]),
		wf_ref<RationalFlux>(g, *this) {
		boost::mt19937 mt; //this is randomness
		boost::uniform_real<> r01(-1, 1); //this is the distribution
		boost::variate_generator<boost::mt19937, boost::uniform_real<> > die_real(mt, r01); //glue randomness with distribution
		double sum2 = 0;
		for(size_t i = 0; i != num_elements(); ++i) {
			origin()[i] = std::complex<double>(die_real(), die_real());
			sum2 = norm(origin()[i]);
		}
		for(size_t i = 0; i != num_elements(); ++i) {
			origin()[i] /= std::sqrt(sum2);
		}
	}
};

template<int IntegerFlux>
struct wf_integer {
	typedef wf<static_rational<IntegerFlux> > type;
};


template<>
field<double> wf_ref<static_rational<0> >::density() const {
	field<double> ret(grid_);
	multi_array<complex, 3> wfxyz(ranges(wfXYZ_ref_));
	plan(wfXYZ_ref_, wfxyz, backward).execute_normalize();
	double const wf_norm = sqrt(ret.data_.num_elements() / ret.domain().volume()); // sqrt(1./((L/grid)^3 )) correct??
	for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
		ret.data_.origin()[i] = norm(wfxyz.origin()[i] / wf_norm);
	}
	return ret;
}

template<>
field<double> wf_ref<static_rational< 1 > >::density() const {
	field<double> ret(grid_);
	const_multi_array_ref<complex, 2> wfXZh(wfXYZ_ref_.origin(), extents[wfXYZ_ref_.shape()[0]][wfXYZ_ref_.shape()[1]*wfXYZ_ref_.shape()[2]]);
	multi_array<complex, 2>       wfXzh(shape(wfXZh));
	plan(-1, wfXZh, wfXzh, backward).execute_normalize();
	multi_array_ref<complex, 3>   wfXYz(wfXzh.origin(), shape(wfXYZ_ref_));
	multi_array<complex, 3>       wfxyz(shape(wfXYZ_ref_));
	plan(+2, wfXYz, wfxyz, backward).execute_normalize();
	double sum = 0;
	double const wf_scale = sqrt(wfxyz.num_elements() / ret.domain().volume()); // sqrt(1./((L/grid)^3 ))
	for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
		ret.data_.origin()[i] = norm(wfxyz.origin()[i] * wf_scale);
		sum += ret.data_.origin()[i];
	}
	//clog<<"sum is "<<sum<<endl;
	//assert(fabs(sum-1.)<0.0000001); //all states are normalized to one in their fft grid/array representation
	return ret;
}

template<>
field<complex> wf_ref<static_rational< 1 > >::wavefunction() const {
	field<complex> ret(grid_);
	const_multi_array_ref<complex, 2> wfXZh(wfXYZ_ref_.origin(), extents[wfXYZ_ref_.shape()[0]][wfXYZ_ref_.shape()[1]*wfXYZ_ref_.shape()[2]]);
	multi_array<complex, 2>       wfXzh(shape(wfXZh));
	plan(-1, wfXZh, wfXzh, backward).execute_normalize();
	multi_array_ref<complex, 3>   wfXYz(wfXzh.origin(), shape(wfXYZ_ref_));
	multi_array<complex, 3>       wfxyz(shape(wfXYZ_ref_));
	plan(+2, wfXYz, wfxyz, backward).execute_normalize();
	double sum = 0;
	double const wf_scale = sqrt(wfxyz.num_elements() / ret.domain().volume()); // sqrt(1./((L/grid)^3 ))
	for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
		ret.data_.origin()[i] = wfxyz.origin()[i] * wf_scale; //could be optimized by putting ret.data_ inside plan(+2,...);
		sum += norm(ret.data_.origin()[i]);
	}
	//CLOG(sum/pow(wf_scale,2));
	assert(fabs(sum / pow(wf_scale, 2) - 1.) < 0.0000001); //all states are normalized to one in their fft grid/array representation
	return ret;
}

using boost::optional;
field<vector> magnetic(
	field <vector> const& current,
	optional<vector> position = optional<vector>(),
	vector* value = 0 
){
	field<vector> bret(current);
	multi_array<array<complex, 3>, 3> jreal(shape(current.data_));
	for(unsigned i = 0; i != bret.data_.num_elements(); ++i) {
		jreal.origin()[i][0] = current.data_.origin()[i][0];
		jreal.origin()[i][1] = current.data_.origin()[i][1];
		jreal.origin()[i][2] = current.data_.origin()[i][2];
	}
	multi_array<array<complex, 3>, 3> jrcpr(shape(jreal));
	fftw3::plan(jreal, jrcpr, fftw3::forward).execute_normalize();
	complex I(0, 1); //make it a class
	double b1, b2, b3;
	if(value) {
		(*value)[0] = (*value)[1] = (*value)[2] = 0.0;
		b1 = current.domain().b(0)[0];
		b2 = current.domain().b(1)[1];
		b3 = current.domain().b(2)[2];
    }	
	boost::wrap_array_view<multi_array_ref<array<complex, 3>, 3> > wjrcpr(jrcpr);
	multi_array< array<complex, 3>, 3> brcpr(shape(jrcpr));
	boost::wrap_array_view<multi_array_ref<array<complex, 3>, 3> > wbrcpr(brcpr);
	for(int i0 = wbrcpr.index_bases()[0]; i0 != (int)(wbrcpr.index_bases()[0] + wbrcpr.shape()[0]); ++i0) {
		for(int i1 = wbrcpr.index_bases()[1]; i1 != (int)(wbrcpr.index_bases()[1] + wbrcpr.shape()[1]); ++i1) {
			for(int i2 = wbrcpr.index_bases()[2]; i2 != (int)(wbrcpr.index_bases()[2] + wbrcpr.shape()[2]); ++i2) {
				//In[181]:= Cross[{gx, gy, gz}, {Jx, Jy, Jz}]
				//Out[181]= {gy Jz - gz Jy, gz Jx - gx Jz, gx Jy - gy Jx}
				//boost::array<int, 3> idx={{i0,i1,i2}};
				double gx, gy, gz, g2;
				gx = i0*b1; gy = i1*b2; gz = i2*b3; g2 = gx*gx + gy*gy + gz*gz;
				wbrcpr[i0][i1][i2][0] = (g2 != 0.) ? (I * (gy * wjrcpr[i0][i1][i2][2] - gz * wjrcpr[i0][i1][i2][1]) / g2) : (0.);
				wbrcpr[i0][i1][i2][1] = (g2 != 0.) ? (I * (gz * wjrcpr[i0][i1][i2][0] - gx * wjrcpr[i0][i1][i2][2]) / g2) : (0.);
				wbrcpr[i0][i1][i2][2] = (g2 != 0.) ? (I * (gx * wjrcpr[i0][i1][i2][1] - gy * wjrcpr[i0][i1][i2][0]) / g2) : (0.);
				if(position){assert(value);
					(*value)[0]+=real(wbrcpr[i0][i1][i2][0]*exp(I*(gx*(*position)[0]+gy*(*position)[1]+gz*(*position)[2])));
					(*value)[1]+=real(wbrcpr[i0][i1][i2][1]*exp(I*(gx*(*position)[0]+gy*(*position)[1]+gz*(*position)[2])));
					(*value)[2]+=real(wbrcpr[i0][i1][i2][2]*exp(I*(gx*(*position)[0]+gy*(*position)[1]+gz*(*position)[2])));
				}
			}
		}
	}
	if(position){assert(value);//normalization
		(*value)[0]/=sqrt(wbrcpr.shape()[0]*wbrcpr.shape()[0]*wbrcpr.shape()[0]);
		(*value)[1]/=sqrt(wbrcpr.shape()[0]*wbrcpr.shape()[0]*wbrcpr.shape()[0]);
		(*value)[2]/=sqrt(wbrcpr.shape()[0]*wbrcpr.shape()[0]*wbrcpr.shape()[0]);
	}
	multi_array< array<complex, 3>, 3> breal(shape(brcpr));
	fftw3::plan(brcpr, breal, fftw3::backward).execute_normalize();
	for(unsigned i = 0; i != bret.data_.num_elements(); ++i) {
		bret.data_.origin()[i][0] = real(breal.origin()[i][0]);
		bret.data_.origin()[i][1] = real(breal.origin()[i][1]);
		bret.data_.origin()[i][2] = real(breal.origin()[i][2]);
	}
	return bret;
}

template<>
field<vector> wf_ref<static_rational<1> >::current() const {
	const_multi_array_ref<complex, 2> wfXZh(
	    wfXYZ_ref_.origin(),
	    extents[shape(wfXYZ_ref_)[0]][shape(wfXYZ_ref_)[1]*shape(wfXYZ_ref_)[2]]
	);
	multi_array<complex, 2> Jz_wfXZh(wfXZh);
	multi_array<complex, 2> Jz_wfXzh(shape(wfXZh));
	multi_array<complex, 2> wfXzh(shape(wfXZh));
	multi_array<complex, 2> Jx_wfXzh(wfXzh);
	multi_array<complex, 2> Jy_wfXzh(wfXzh);

	bool wrap[] = {true, true};
	boost::wrap_array_view<multi_array_ref<complex, 2> > WJz_wfXZh(Jz_wfXZh, wrap);
	{
		for(int i = WJz_wfXZh.index_bases()[0]; i != (int)(WJz_wfXZh.index_bases()[0] + WJz_wfXZh.shape()[0]); ++i) {
			for(int j = WJz_wfXZh.index_bases()[1]; j != (int)(WJz_wfXZh.index_bases()[1] + WJz_wfXZh.shape()[1]); ++j) {
#ifndef _SHORTFFT	
				WJz_wfXZh[i][j] *= (j / (double)shape(wfXYZ_ref_)[1]) * domain().b(2)[2];
#else
				WJz_wfXZh[i][j] *= (j / (double)shape(wfXYZ_ref_)[1]) * domain().b(2)[2] * 2.0;
#endif
			}
		}
	}

#ifndef _SHORTFFT
	plan(-1, Jz_wfXZh, Jz_wfXzh, fftw3::backward).execute_normalize();
	plan(-1, wfXZh, wfXzh, fftw3::backward).execute_normalize();
#else
	FourierTransform ft( (boost::array<unsigned,3>) 
			{{shape(wfXYZ_ref_)[0],shape(wfXYZ_ref_)[1],shape(wfXYZ_ref_)[2]}}
			);
	ft.z_backward_mpbc(wfXZh, wfXzh);
	ft.z_backward_mpbc(Jz_wfXZh, Jz_wfXzh);
#endif

	bool wrap2[] = {true, false};
	{
		boost::wrap_array_view<multi_array_ref<complex, 2> > WJx_wfXzh(Jx_wfXzh, wrap2);
		for(int i = WJx_wfXzh.index_bases()[0]; i != (int)(WJx_wfXzh.index_bases()[0] + WJx_wfXzh.shape()[0]); ++i) {
			for(int j = WJx_wfXzh.index_bases()[1]; j != (int)(WJx_wfXzh.index_bases()[1] + WJx_wfXzh.shape()[1]); ++j) {
				WJx_wfXzh[i][j] *= (i / (double)WJx_wfXzh.shape()[0]) * domain().b(0)[0];
			}
		}
	}
	bool wrap3[] = {true, true};
	boost::wrap_array_view<multi_array_ref<complex, 2> > WJy_wfXzh(Jy_wfXzh, wrap3);
	{
		for(int i = WJy_wfXzh.index_bases()[0]; i != (int)(WJy_wfXzh.index_bases()[0] + WJy_wfXzh.shape()[0]); ++i) {
			for(int j = WJy_wfXzh.index_bases()[1]; j != (int)(WJy_wfXzh.index_bases()[1] + WJy_wfXzh.shape()[1]); ++j) {
				WJy_wfXzh[i][j] *= (j / (double)shape(wfXYZ_ref_)[2]) * domain().b(1)[1];
			}
		}
	}	

	multi_array_ref<complex, 3> Jx_wfXYz(Jx_wfXzh.origin(), shape(wfXYZ_ref_));
	multi_array<complex, 3> Jx_wfxyz(shape(wfXYZ_ref_));
	plan(+2, Jx_wfXYz, Jx_wfxyz, fftw3::backward).execute_normalize();
	multi_array_ref<complex, 3> wfXYz(wfXzh.origin(), shape(wfXYZ_ref_));
	multi_array_ref<complex, 3> Jy_wfXYz(Jy_wfXzh.origin(), shape(wfXYz));
	multi_array<complex, 3> Jy_wfxyz(shape(wfXYZ_ref_));
	plan(+2, Jy_wfXYz, Jy_wfxyz, fftw3::backward).execute_normalize();
	multi_array_ref<complex, 3> Jz_wfXYz(Jz_wfXzh.origin(), shape(wfXYZ_ref_));
	multi_array<complex, 3> Jz_wfxyz(shape(wfXYZ_ref_));
	plan(+2, Jz_wfXYz, Jz_wfxyz, fftw3::backward).execute_normalize();
	boost::array<multi_array<double, 3> , 3> Jreal;
	for(unsigned i = 0; i != Jreal.size(); ++i) {
		Jreal[i].resize(shape(wfXYZ_ref_));
	}
	multi_array    <complex, 3> wfxyz(shape(wfXYZ_ref_));
	plan(+2, wfXYz, wfxyz, fftw3::backward).execute_normalize();
	for(unsigned i = 0; i < wfxyz.num_elements(); ++i) {
		Jreal[0].origin()[i] = real(conj(wfxyz.origin()[i]) * Jx_wfxyz.origin()[i]);
		Jreal[1].origin()[i] = real(conj(wfxyz.origin()[i]) * Jy_wfxyz.origin()[i]);
		Jreal[2].origin()[i] = real(conj(wfxyz.origin()[i]) * Jz_wfxyz.origin()[i]);
	}
	field<vector> J(grid_);
#ifndef _SHORTFFT
	double const wf_scale2 = (wfxyz.num_elements() / J.domain().volume()); // sqrt(1./((L/grid)^3 ))
#else
	/* need to double check the rationale of this normalization factor */
	double const wf_scale2 = (wfxyz.num_elements() * 2 / J.domain().volume());
#endif
	for(unsigned i = 0; i < wfxyz.num_elements(); ++i) {
		J.data_.origin()[i][0] = Jreal[0].origin()[i] * wf_scale2;
		J.data_.origin()[i][1] = Jreal[1].origin()[i] * wf_scale2;
		J.data_.origin()[i][2] = Jreal[2].origin()[i] * wf_scale2;
	}
	return J;
}

template<>
field <vector> wf_ref<static_rational<0> >::current() const {
	field <vector> Jreal(grid_);
	return Jreal;
}

field<vector> current(
    Wavefunction const& Wf_,
    optional<unsigned> occupation = optional<unsigned>()
) {
#ifndef _REGULARGDISTRIB
	assert(0); //grid must be regular G
	//BOOST_STATIC_ASSERT(0);
#endif
	field <
	//array<double, 3>
	vector
	> ret(get_grid(Wf_));
	for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
		ret.data_.origin()[i][0] = 0;
		ret.data_.origin()[i][1] = 0;
		ret.data_.origin()[i][2] = 0;
	}
	assert(ret.data_.origin()[20][2] == 0.);
	for(unsigned istate = 0; istate != (unsigned)Wf_.sd(0, 0)->nstloc(); ++istate) {
		double used_occ = occupation ? (istate < *occupation ? 1. : 0.) : Wf_.sd(0, 0)->occ()[istate];
		std::clog << "current for state " << istate << " using occ " << used_occ << std::endl;
		field<vector> current = qbox::wf_ref<static_rational<
#ifndef _MPBC
			0
#else
			1
#endif
		> > (Wf_, istate).current();
		for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
			typedef 
				linear_space<vector, std::plus<vector>,	double,	stdx::multiplies<double, vector, vector> > 
				space;
			space cdoi = current.data_.origin()[i];
			space cdoi_w = (used_occ) * cdoi;
			(space&)ret.data_.origin()[i]
			+= cdoi_w;
		}
	}
	return ret;
}

#ifndef _MPBC
#define FLUX 0
#else
#define FLUX 1
#endif

field<double> density(qbox::wf<static_rational<FLUX> > const& wf);

field<double> density(qbox::wf<static_rational<FLUX> > const& wf){
	return wf.density();
}

field<complex> wavefunction(
	Wavefunction const& Wf, 
	unsigned idx_state
){
	return qbox::wf_ref<static_rational< FLUX > >(Wf, idx_state, 0, 0).wavefunction();
}

using boost::optional;
field<double> density(
    Wavefunction const& Wf_,
    optional<unsigned> occupation = optional<unsigned>()
) {
#ifndef _REGULARGDISTRIB
	assert(0); //grid must be regular g to use this function
	//BOOST_STATIC_ASSERT(0);
#endif
	field<double> ret(get_grid(Wf_)); /* compute charge */
	for(unsigned ispin = 0; ispin < numeric_cast<unsigned>(Wf_.nspin()) ; ispin++) {
		for(unsigned ikp = 0; ikp < numeric_cast<unsigned>(Wf_.nkp()); ikp++) {
			for(unsigned istate = 0; istate != numeric_cast<unsigned>(Wf_.sd(ispin, ikp)->nstloc()); ++istate) {
				double used_occ = occupation ? (istate < *occupation ? 1. : 0.) : Wf_.sd(ispin, ikp)->occ()[istate];
				std::clog << "density for state " << istate << " using occ " << used_occ << std::endl;
				field<double> density = qbox::wf_ref<static_rational< FLUX > >(Wf_, istate, ispin, ikp).density();
				for(unsigned i = 0; i != ret.data_.num_elements(); ++i) {
					ret.data_.origin()[i] += used_occ * Wf_.weight(ikp) * density.data_.origin()[i];
				}
			}
		}
	}
	return ret;
}

qbox::temperature kinetic_temperature(AtomSet const& atoms) {
	std::vector<std::vector<double> > v;
	atoms.get_velocities(v);
	qbox::energy ekin = 0.* atomic::E_h;
	unsigned ndofs = 0;
	//pmass_[is] = atoms_.species_list[is]->mass() * 1822.89;
	quantity<atomic::mass> m_u = 1822.89 * atomic::m_e; /* = dalton*/
	for(unsigned is = 0; is < v.size(); is++) {
		for(unsigned i = 0; i < v[is].size(); i++) {  //warning atom 1 vx vy vz, atom 2 vx vy vz
			quantity<atomic::velocity> v_is_i = v[is][i] * (atomic::a_0 / atomic::time_unit);
			qbox::energy dekin =
			    0.5 *
			    atoms.species_list[is]->mass() /*atomic-mass*/ * m_u
			    * v_is_i * v_is_i;
			ekin = ekin + dekin;
			ndofs++;
		}
	}
	qbox::temperature ret((2. * ekin / (double)ndofs) / si::constants::codata::k_B.value());
	return ret;
}
session const& session::export_xyz(boost::filesystem::path const& p, std::string comment /*="# -1"*/) const {
	history_ << "# atoms exported as xyz file '" << p << "'" << endl;
	std::ostringstream oss;
	if(comment=="# -1") comment="# snapshot T_ion =" + boost::lexical_cast<std::string>(qbox::kinetic_temperature(this->atoms()));
	this->export_xyz(oss, comment);
	boost::filesystem::ofstream ofs(p);
	ofs << oss.str() << std::flush;
	return *this;
}
session& session::import_xyz(boost::filesystem::path const& p, std::map<std::string, boost::filesystem::path> const& species){
	boost::filesystem::ifstream ifs(p);
	xyz::frame const f = xyz::frame::load(p);
	double atob(quantity<atomic::length>(1.*nonsi::angstrom)/atomic::bohr);
	if(f){
		cell c(
			   vector((*f)[0][0]*atob, (*f)[0][1]*atob, (*f)[0][2]*atob), 
			   vector((*f)[1][0]*atob, (*f)[1][1]*atob, (*f)[1][2]*atob), 
			   vector((*f)[2][0]*atob, (*f)[2][1]*atob, (*f)[2][2]*atob)
		);
		this->set_cell(c);
		//this->set_refcell(c*1.05);
	}
	for(std::map<std::string, boost::filesystem::path>::const_iterator it=species.begin(); it!=species.end(); ++it){
		this->species(it->first, it->second);
	};
	for(unsigned i=0; i!=f.size(); ++i){
		if(species.find(f[i])==species.end()) throw std::runtime_error("don't know how to translate species with name "+f[i]+" to qbox");
		std::string name = f[i]+ boost::to_string(i+1);
		qbox::vector v(
					   f[i][0].value()*atob, 
					   f[i][1].value()*atob, 
					   f[i][2].value()*atob
					   );
		this->atom(name, (std::string const&)f[i], v);
	}
	return *this;
}
//namespace qbox
}
#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 cindent: */
