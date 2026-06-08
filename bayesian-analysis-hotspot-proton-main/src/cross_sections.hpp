#ifndef cross_sections_hpp
#define cross_sections_hpp

class	xsections{
	public:
		/* 
			(--) code_Nh = 0: without hotspot number fluctuation 
			(--) code_Nh = 1: with hotspot number fluctuation. 
		*/
		double dsigma_coherent(int code_Nh, double Q2, double W, double t, model_params mparams);
		//double dsigma_coherent_rcut(int code_Nh, double Q2, double W, double t, model_params mparams);
		double dsigma_ci(int code_Nh, double Q2, double W, double t, model_params mparams);
		double dsigma_ci_rcut(int code_Nh, double Q2, double W, double t, model_params mparams);
		double dsigma_ci_rcut_red(int code_Nh, double Q2, double W, double t, model_params mparams);
};


#endif /*cross-sections*/
