#include "metis1.h"




int metis11() {
	idx_t nn, ne, num_node;
	int  num_element, num_surface;
	int i, j, cc;
	string filename = "BulkInfoooo";
	string filename1 = ".mphtxt";
	string result1 = filename + "_epart.txt";
	extern Eigen::MatrixXd zb;
	extern Eigen::MatrixXi ver;
	Eigen::VectorXi element;
	ifstream infile;
	ofstream ofs;
	ofs.precision(16);
	infile.open(filename + filename1);
	//infile >> num_domain;
	//num_domain = ii;
	cout << "num_domain=" << num_domain << endl;
	infile >> num_node;
	cout << "num_node=" << num_node << endl;
	/***************************************************************************/
	/*                                                                         */
	/*             read the mesh file from the hard disk	                   */
	/*                                                                         */
	/***************************************************************************/
	zb.resize(num_node, 3);
	for (cc = 0; cc < num_node; cc++) {
		for (j = 0; j < 3; j++) {
			infile >> zb(cc, j);
		}
	}
	infile >> num_element;
	ne = num_element;
	cout << "num_element=" << num_element << endl;
	ver.resize(num_element, 12);
	for (int el = 0; el < num_element; el++) {
		for (int i = 0; i < 4; i++) {
			infile >> ver(el, i);
			ver(el, i)++;
		}
	}
	element.resize(num_element);
	for (int i = 0; i < num_element; i++) {
		infile >> element[i];
	}
	infile.close();
	nn = 4 * ne;
	nn = num_node;
	idx_t ncommon = 3;
	idx_t nParts = num_domain;
	idx_t objval;
	vector<idx_t> epart(ne, 0);
	vector<idx_t> npart(nn, 0);
	vector<idx_t> count(nParts, 0);
	vector<idx_t> eptr(0);
	vector<idx_t> eind(0);
	idx_t options1[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options1);
	options1[METIS_OPTION_NUMBERING] = 0;
	for (i = 0; i < ne; i++) {
		eptr.push_back(4 * i);
		for (j = 0; j < 4; j++) {
			eind.push_back(ver(i, j) - 1);
		}
	}
	eptr.push_back(4 * ne);
	int ret = METIS_PartMeshDual(&ne, &nn, eptr.data(), eind.data(), NULL, NULL, &ncommon, &nParts, NULL, NULL, &objval, epart.data(), npart.data());
	for (i = 0; i < ne; i++) {
		count[epart[i]]++;
	}
	ofs.open(result1, ios::out);
	for (i = 0; i < ne; i++) {
		ofs << epart[i] + 1 << '\n';
	}
	return 0;
}