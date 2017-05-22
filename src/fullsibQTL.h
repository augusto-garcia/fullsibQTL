void reorg_genoprob(int n_ind, int n_pos, int n_gen, 
		    double *genoprob, double ****Genoprob);

void allocate_double(int n, double **vector);

void allocate_dmatrix(int n_row, int n_col, double ***matrix);

void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod);

double log_add(double x, double y);

void matmult(double *result, double *a, int nrowa,
             int ncola, double *b, int ncolb);

void R_scan_qtl(int *n_ind, int *n_pos, int *n_gen, int *colin_tp, 
		double *genoprob, double *addcov, int *nc_Addcov,
		double *gamma, double *addcov_ls, double *pheno, 
		double *result, int *maxit, double *tol, int *verbose);


void scan_qtl(int n_ind, int n_pos, int *n_gen, double ***Genoprob,
	      double *pheno,  double **Addcov, int nc_Addcov, 
	      double *gamma, double **Addcov_ls, double *result, 
	      int maxit, double tol, int *colin_tp, int verbose);

double em(int n_ind, int n_gen, double **Geno1pos, double *pheno, 
	  double **Addcov, int nc_Addcov, double *gamma,  double **Addcov_ls, 
	  int maxit, double tol, double *effects,
	  int hypot, int colin_tp, int verbose);

void contrast_used (int colin_tp, double **contr);

void restr_effects(int n_gen, int hypot, double *effects);

int H0test (int colin_tp);



void R_char_qtl(int *n_ind, int *n_gen, int *colin_tp, 
		double *genoprob, double *addcov, int *nc_addcov, 
		double *gamma, double *addcov_ls, double *pheno,
		double *result, int *maxit, double *tol, int *verbose);

void char_qtl(int n_ind, int n_gen, double **Genoprob,
	      double *pheno,  double **Addcov, int nc_Addcov, 
	      double *gamma, double **Addcov_ls, double *result, 
	      int maxit, double tol, int colin_tp, int verbose);


void qtl_3effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose);

void qtl_2effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose);

void qtl_1effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose);
