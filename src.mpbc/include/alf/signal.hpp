#include <csignal>
#include <iostream>
//#include<boost/signal.hpp>

namespace system_signal{
	namespace terminal_interrupt{
		static bool flag_=false;
		void handle_helper(int sig){
			std::clog
				<<"catched signal: "
				<<sig<<" ... flag set and will be handled, hit again to confirm and quit immediately!"
				<<"\n"<<sig;
			(void)signal(SIGINT, SIG_DFL);
			flag_ = true;
		}
		void capture(){
			std::clog<<"control-C handled"<<std::endl;
			(void)signal(SIGINT,handle_helper);
		}
		void release(){
			std::clog<<"control-C not handled"<<std::endl;
			(void)signal(SIGINT,SIG_DFL);
		}
		bool flag(){
			return flag_;
		}
	}
}
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:

