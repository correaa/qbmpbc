#ifndef GSL_ERROR_HPP__
#define GSL_ERROR_HPP__
#include<gsl/gsl_errno.h>
#include<string>
#include<boost/lexical_cast.hpp>
#include<stdexcept>
#include<cassert>
#include<iostream> //debug
namespace gsl{
namespace error{
//	enum core;
	enum code{
		input_domain_error = GSL_EDOM,     // 1
		output_range_error = GSL_ERANGE,   // 2
		no_memory          = GSL_ENOMEM,   // 8
		invalid_argument   = GSL_EINVAL,   // 4
		bad_function       = GSL_EBADFUNC, // 9
		failure            = GSL_FAILURE   // -1
	};
	struct status{
		code c_;
		status(code const& c) : c_(c){}
		status(int const& c) : c_((code)c){}
		operator code() const{return c_;}
		operator std::string() const{return gsl_strerror(c_);}
	};
	struct exception : std::runtime_error{
		int const code_;
		exception(int gsl_errno, std::string const& c="no detailed description") : std::runtime_error(std::string(gsl_strerror(gsl_errno)) + ": " + c), code_(gsl_errno){
			if(gsl_errno==0) std::clog << "error code 0 is not an error!" << std::endl;
		}
		int code() const{return code_;}
	};
	void handler(const char * reason, const char * file, int line, int gsl_errno){
		throw exception(
				gsl_errno, 
			  std::string(reason)
			+ " at " + std::string(file)+":"+boost::lexical_cast<std::string>(line) 
			+ " error code "+boost::lexical_cast<std::string>(gsl_errno) + ": " + std::string(gsl_strerror(gsl_errno))
		);
	}
	gsl_error_handler_t* set_handler(void (*a)(const char*, const char*, int, int)){
		assert(a);
		return gsl_set_error_handler(a);
	}
	gsl_error_handler_t* unset_handler(){
		return gsl_set_error_handler_off();
	}
	static gsl_error_handler_t* const native = set_handler(&handler);
}}
#endif

