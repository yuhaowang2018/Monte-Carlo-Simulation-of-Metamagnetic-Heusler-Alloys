#include "lattice.h"
#include "MersenneTwister.h"



lattice::lattice()
{
    REPNOS = 0;
    ORGCELLREPLICA = 6;
    c2aRatio = 1.0;
    STRUCTCUBDIS = 1.0;
    STRUCTTETDIS = 1.0;
    STRDATATYPE = 2;
    Kaniscub =0;
    Kanistet =0;
    Structure   =   "L21";
    DIRAniso    =   "ZDIR";
    POLY=0;

}


lattice::~lattice()
{
    //dtor
}

void lattice::setatoms(vector<atom> atoms)
{
  _atoms=atoms;

}

vector<atom> lattice::getatoms()
{   return _atoms;

}

void lattice::AtomSubstitute(int Atin, int Atout, int COUNTREP)
{
    int COUNTING = 0, ATOMSEL;
    while ( COUNTING != COUNTREP )
    {
        ATOMSEL = RNDMNO.randInt( _SIMDOMAINATOMS - 1 );
        if ( _SIMDOMAINDATA[ATOMSEL].getnumber() == Atout )
        {   atom aa=_SIMDOMAINDATA[ATOMSEL];
            double x=_SIMDOMAINDATA[ATOMSEL].getx();
            double y=_SIMDOMAINDATA[ATOMSEL].gety();
            double z=_SIMDOMAINDATA[ATOMSEL].getz();
            _SIMDOMAINDATA[ATOMSEL]=_atoms[Atin];
            _SIMDOMAINDATA[ATOMSEL].setcoordinate(x,y,z);
            atom bb=_SIMDOMAINDATA[ATOMSEL];
            COUNTING ++;
        }
    }
}

void lattice::CreateSimDomain(modelstats modelnum)
{ fstream sdFILE;
   int REPNOS, initAtoms, Node;
    double  **initAtCood;
    /**Read SimDomain Data from the file***/
    if ( ! ( strcmp (modelnum.DomainFileFlag.c_str(), "READ" ) ) )
    {

    }
    else
    {
        sdFILE.open ( modelnum.DomainFileName.c_str(), ios::out );
        if ( !sdFILE.is_open () )
        {
            cout<<"Unable to open "<<modelnum.DomainFileName<<" to write ... Exiting\n";
            exit(0);
        }
    /*** If data is not read from a file ***/
         if (!(strcmp(Structure.c_str(), "B2"))) {
        REPNOS = EAATOMTOMS-2;
        initAtoms = 2;
        //funit = 2;
        initAtCood = new double *[initAtoms];
        for (int i=0; i<initAtoms; i++) initAtCood[i] = new double [4];
        double B2Cood[2][4] = {{0, 0.0, 0.0, 0.0}, {1, 0.5, 0.5, 0.5}};
        for (int i=0; i<initAtoms; i++) for (int j=0; j<4; j++) initAtCood[i][j] = B2Cood[i][j];
    }
    else if (!(strcmp(Structure.c_str(), "L21")))   {
        REPNOS = EAATOMTOMS-3;
        initAtoms = 16;
        //funit = 4;
        initAtCood = new double *[initAtoms];
        for (int i=0; i<initAtoms; i++) initAtCood[i] = new double [4];
        double L21Cood[16][4] = {{2, 0.0, 0.0, 0.0}, {2, 0.0, 0.5, 0.5}, {2, 0.5, 0.5, 0.0}, {2, 0.5, 0.0, 0.5},
                                 {1, 0.0, 0.5, 0.0}, {1, 0.0, 0.0, 0.5}, {1, 0.5, 0.0, 0.0}, {1, 0.5, 0.5, 0.5},
                                 {0, 0.25, 0.25, 0.25}, {0, 0.25, 0.75, 0.25}, {0, 0.25, 0.25, 0.75}, {0, 0.25, 0.75, 0.75},
                                 {0, 0.75, 0.25, 0.25}, {0, 0.75, 0.75, 0.25}, {0, 0.75, 0.25, 0.75}, {0, 0.75, 0.75, 0.75}};
        for (int i=0; i<initAtoms; i++) for (int j=0; j<4; j++) initAtCood[i][j] = L21Cood[i][j];
    }

    _SIMDOMAINATOMS = initAtoms*ORGCELLREPLICA*ORGCELLREPLICA*ORGCELLREPLICA;
    _SIMDOMAINDATA  =new atom [_SIMDOMAINATOMS];
    for (int z=0; z<ORGCELLREPLICA; z++) for (int y=0; y<ORGCELLREPLICA; y++) for (int x=0; x<ORGCELLREPLICA; x++) for (int eaAtom = 0; eaAtom<initAtoms; eaAtom++)  {
        Node = eaAtom+(x+y*ORGCELLREPLICA+z*ORGCELLREPLICA*ORGCELLREPLICA)*initAtoms;
        int number=initAtCood[eaAtom][0];
        atom aa=_atoms[number];
        _SIMDOMAINDATA[Node]=aa;

        _SIMDOMAINDATA[Node].setcoordinate(initAtCood[eaAtom][1] + x,initAtCood[eaAtom][2] + y,initAtCood[eaAtom][3] + z);
        atom bb=_SIMDOMAINDATA[Node];

    }

        /** Count the number of atoms **/



        for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
        {
            int count1=_atoms[int(_SIMDOMAINDATA[i].getnumber())].getcount();
            count1++;
            _atoms[int(_SIMDOMAINDATA[i].getnumber())].setcount(count1);
        }
        cout<<"\nATOM COUNT BEFORE REPLACEMENT OF ATOMS ... \n";
        for ( int i = 0; i < EAATOMTOMS; i++)
        {
            cout<<_atoms[i].name<<"\t"<<_atoms[i].getcount()<<endl;
        }

        /** Replace Atoms **/
        if ( REPNOS > 0 )
        {
            int ATOM1, ATOM2, COUNTREP;
            for ( int i = 0; i < REPNOS; i++ )
            {
                ATOM1 = REPDATA[i][0];
                ATOM2 = REPDATA[i][1];
                COUNTREP = int ( REPDATA[i][2] * _atoms[ATOM2].getcount() + 0.5 );  // 0.5 to get nearest integer
                cout << "\nSubstuting "<<COUNTREP<<" ATOMS of "<<_atoms[ATOM2].name<<" with "<<_atoms[ATOM1].name<<" ..... "<<endl;
                AtomSubstitute ( ATOM1, ATOM2, COUNTREP );
            }
        }

        for( int i = 0; i < _SIMDOMAINATOMS; i++ )
        {
           _atoms[int(_SIMDOMAINDATA[i].getnumber())].setcount(0);

        }

        for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
        {
            int count1=_atoms[int(_SIMDOMAINDATA[i].getnumber())].getcount();
            count1++;
            _atoms[int(_SIMDOMAINDATA[i].getnumber())].setcount(count1);
        }

        cout<<"\nATOM COUNT after replace SIM DOMAIN : \n";
        for ( int i = 0; i < EAATOMTOMS; i++)
        {
            cout<<_atoms[i].name<<"\t"<<_atoms[i].getcount()<<endl;
        }

        /*** Write sim domain data to file ***/

        sdFILE<<"TotalAtoms \t"<<_SIMDOMAINATOMS<<endl;
        for ( int i = 0; i < EAATOMTOMS; i++)
        {
            sdFILE<<_atoms[i].name<<"\t"<<_atoms[i].getcount()<<endl;
        }
        for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
        {
            sdFILE<<_SIMDOMAINDATA[i].name<<"\t"<<_SIMDOMAINDATA[i].getx()<<"\t"<<_SIMDOMAINDATA[i].gety()<<"\t"<<_SIMDOMAINDATA[i].getz();
            sdFILE<<"\n";
        }

        /******************************
        ERASE Memory Alloc not required
        ******************************/
//        for ( int i = 0; i < ATOMSINCELL; i++ )
//        {
//            delete [] ATCOORDATA[i];
//        }
//        for ( int i = 0; i < ORGATOMNO; i++ )
//        {
//            delete [] ORGATOMDATA[i];
//        }
//        delete [] ORGATOMDATA;
//        delete [] ATCOORDATA;
    }
    sdFILE.close();
}




void lattice::CreateMagNNList (string NAME, int INITJ, double **DATAINITJ, double _c2aRatio)
{
    /****************************************
    MAGCUBICPOINTS -> INITJ - External
        Intial J Parameter between atoms
    TEMPMAGCUBICPTS -> FINALJ
        Final J points after adding symmetries
    MAGCUBICDATA -> DATAINITJ - External
        Data list of atoms, distance and J values
    TEMPMAGCUBICDATA -> DATAFINALJ
        Data list of atoms, distance and J values
        after adding symmetries
    SIMDOMAINATOMS -> SIMDOMAINATOMS - External
        Total atoms in the simualtion domain
    SIMDOMAINDATA -> SIMDOMAINDATA - External
        Atomic co-ordinate  data of all the atoms
        of simualtion domain
    MAGCUBCOUNT -> COUNTER
        Counter to generate final J data
    ORGCELLREPLICA -> ORGCELLREPLICA - External
        Repetations of the original cell in all directions
    MAGCUBICNNLIST -> NNLIST - External
        Neighbor list generated by the function
    ****************************************/

    /** Check magnetic data **/
    int ATOMi, ATOMj, FINALJ;
    FINALJ = INITJ;
    for ( int i = 0 ; i < INITJ; i++ )
    {
        ATOMi = int(DATAINITJ[i][0]);
        ATOMj = DATAINITJ[i][1];
        if ( ATOMi != ATOMj )
        {
            for ( int j = 0 ; j < INITJ; j++ )
            {
                if ( DATAINITJ[i][2] == DATAINITJ[j][2] )
                {
                    if ( ATOMj == DATAINITJ[j][0] && ATOMi == DATAINITJ[j][1] )
                    {
                        break;
                    }
                }
                if ( j == INITJ - 1 )
                {
                    FINALJ += 1;
                }
            }
        }
    }
    cout<<"MAGCUBICDATA"<<DATAINITJ<<endl;
    cout<<"List of Chains for "<<NAME<<" neighbor list : "<<FINALJ<<endl;

    double **DATAFINALJ;
    int *COUNTJPERAT;
    COUNTJPERAT = new int [EAATOMTOMS];
    if ( !COUNTJPERAT )
	{
		cout<<"Error allocating memory to COUNTJPERAT in CreateCubicNNLIST ! Exiting ... "<<endl;
		exit(0);
	}
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        COUNTJPERAT[i] = 0;
    }
    DATAFINALJ = new double *[FINALJ];
    for ( int i = 0; i < FINALJ; i++ )
    {
        DATAFINALJ[i] = new double [4];
    }
    if ( !DATAFINALJ )
	{
		cout<<"Error allocating memory to DATAFINALJ in CreateCubicNNLIST ! Exiting ... "<<endl;
		exit(0);
	}

    /** Create new magnetic Cubic parameter lists **/
    int COUNTER = 0;
    for ( int i = 0 ; i < INITJ; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            DATAFINALJ[COUNTER][j] = DATAINITJ[i][j];
        }
        COUNTJPERAT[int(DATAINITJ[i][0])] += 1;
        COUNTER++;
    }
    for ( int i = 0 ; i < INITJ; i++ )
    {
        ATOMi = DATAINITJ[i][0];
        ATOMj = DATAINITJ[i][1];
        if ( ATOMi != ATOMj )
        {
            for ( int j = 0 ; j < INITJ; j++ )
            {
                if ( DATAINITJ[i][2] == DATAINITJ[j][2] )
                {
                    if ( ATOMj == DATAINITJ[j][0] && ATOMi == DATAINITJ[j][1] )
                    {
                        break;
                    }
                }
                if ( j == INITJ - 1 )
                {
                    DATAFINALJ[COUNTER][0] = ATOMj;
                    DATAFINALJ[COUNTER][1] = ATOMi;
                    DATAFINALJ[COUNTER][2] = DATAINITJ[i][2];
                    DATAFINALJ[COUNTER][3] = DATAINITJ[i][3];
                    COUNTJPERAT[ATOMj] += 1;
                    COUNTER += 1;
                }
            }
        }
    }

    cout <<NAME << " J List Before Sorting -> \n";
    for ( int i = 0; i < FINALJ; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            cout << DATAFINALJ[i][j]<<"\t";
        }
        cout<<endl;
    }

    /** Sort the cubic exchange list by
    First -> ATOMi & Second -> Distance **/
    int ATOMi1, ATOMi2, ATOMj1, ATOMj2;
    double DISi, DISj;
    for ( int i = 0; i < FINALJ; i++ )
    {
        for ( int j = i; j < FINALJ; j++ )
        {
            ATOMi1 = DATAFINALJ[i][0];
            ATOMi2 = DATAFINALJ[i][1];
            DISi = DATAFINALJ[i][2];
            ATOMj1 = DATAFINALJ[j][0];
            ATOMj2 = DATAFINALJ[j][1];
            DISj = DATAFINALJ[j][2];
            if ( ATOMj1 < ATOMi1 )
            {
                SwapOneDMatrix ( i, j, 4, DATAFINALJ );
            }
            else if ( ATOMi1 == ATOMj1 )
            {
                if ( ATOMj2 < ATOMi2)
                {
	            SwapOneDMatrix ( i, j, 4, DATAFINALJ );
                }
                else if ( ATOMi2 == ATOMj2 )
                {
                    if ( DISj < DISi )
                    {
                        SwapOneDMatrix ( i, j, 4, DATAFINALJ );
                    }
                }
            }
        }
    }

    cout <<"\n"<<NAME<< " J Exchange List after sorting-> \n";
    for ( int i = 0; i < FINALJ; i++ )
    {
        for(int j=0; j<2;j++)
        {
          cout<<_atoms[int(DATAFINALJ[i][j])].name<<" "<<DATAFINALJ[i][j]<<"\t";
        }
        for ( int j = 2; j < 4; j++ )
        {
            cout << DATAFINALJ[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<< " Count J for each Atom : ";
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        cout << COUNTJPERAT[i] << "\t";
    }
    cout<<"\n";

    /*** Create Cubic Neighbor Lists **/
    cout << "\nCreating "<<NAME<< " Neighbour List ... ";
    cout<<"\n";
    int JCOUNT, JNO;
    double XCOOD, YCOOD, ZCOOD, DISij, ERR, CX, CY, CZ, XHALF, YHALF, ZHALF;
    cout<<ORGCELLREPLICA<<endl;
    XHALF = ORGCELLREPLICA*0.5;
    YHALF = XHALF;
    ZHALF = XHALF*_c2aRatio;
    //NNLIST.resize(_SIMDOMAINATOMS);
    for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
    {
        ATOMi1 = int(_SIMDOMAINDATA[i].getnumber());
        map <double, vector<int>> NNMap;


        JCOUNT = COUNTJPERAT[ATOMi1];
        if ( JCOUNT > 0 )
        {
            /** Code to add J parameter at start **/
            for ( int VAR = 0; VAR < FINALJ; VAR++ )
            {
                if ( ATOMi1 == DATAFINALJ[VAR][0] )
                {

                 for ( int k = 0; k < JCOUNT; k++ )
                    {      vector<int> aa(0);
                        NNMap[DATAFINALJ[VAR+k][3]]=aa;
                    }
                    break;

                 }

            }
          for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
            {
                for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
                {
                    for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                    {
                        XCOOD = _SIMDOMAINDATA[i].getx();
                        YCOOD = _SIMDOMAINDATA[i].gety();
                        ZCOOD = _SIMDOMAINDATA[i].getz()*_c2aRatio;
                        if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                        else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                        if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                        else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                        if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*_c2aRatio;
                        else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*_c2aRatio;
                        for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
                        {
                            if ( j != i )
                            {
                                ATOMj1 = int(_SIMDOMAINDATA[j].getnumber());

                                CX = XCOOD-_SIMDOMAINDATA[j].getx();
                                CY = YCOOD-_SIMDOMAINDATA[j].gety();
                                CZ = ZCOOD-_SIMDOMAINDATA[j].getz()*_c2aRatio;
                                DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));

                                if ( DISij < MAXNNDIS )
                                {   //cout<<"DISij "<<DISij<<" ";
                                    JNO = -1;
                                    //cout<<"FINALJ "<<FINALJ<<endl;
                                    for ( int VAR = 0;  VAR < FINALJ; VAR++ )
                                    {
                                        ATOMi2 = int(DATAFINALJ[VAR][0]);
                                        if ( ATOMi1 == ATOMi2 )
                                        {
                                            JNO += 1;
                                            ATOMj2 = int(DATAFINALJ[VAR][1]);
                                            //cout<<"JNO"<<JNO<<endl;
                                            if ( ATOMj1 == ATOMj2 )
                                            {  //cout<<"ATOMj1 "<<ATOMj1<<endl;

                                                ERR = std::abs(DISij-DATAFINALJ[VAR][2]);

                                                if ( ERR < TOL )
                                                {   //cout<<DATAFINALJ[VAR][3]<<endl;
                                                    NNMap[DATAFINALJ[VAR][3]].push_back(j);
                                                    //cout<<JNO<<endl;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ERR=0;
                }
            }

        }
        if(_c2aRatio==1)
        {
           _SIMDOMAINDATA[i].SetMCNN(NNMap);

        }
        else{
          _SIMDOMAINDATA[i].SetMTNN(NNMap);
        }
    }

    int NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
        {
            if ( _SIMDOMAINDATA[j].getnumber() == i )
            {
                if(_c2aRatio==1){
                cout<<"For Atom \""<<_SIMDOMAINDATA[j].name<<"\" Neighbours number are "<<endl;
                if(COUNTJPERAT[i]==0)
                {cout<<"None"<<endl;}
                map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getMCNN();
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {cout << it->first << " => " << it->second.size() << '\n';}

                NoCount++;
                if (NoCount > 1 ) break;
                }
                else{
                cout<<"For Atom \""<<_SIMDOMAINDATA[j].name<<"\" Neighbours number are "<<endl;
                if(COUNTJPERAT[i]==0)
                {cout<<"None"<<endl;}
                map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getMTNN();
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {std::cout << it->first << " => " << it->second.size()<<'\n';}

                NoCount++;
                if (NoCount > 1 ) break;
                }
            }
        }
    }

    for ( int i = 0; i < FINALJ; i++ )
    {
        delete [] DATAFINALJ[i];
    }
    delete [] DATAFINALJ;
    delete [] COUNTJPERAT;

}



void lattice::MAXStrNNList(string NAME, double _c2aRatio)
{ /***********************************************************




    SIMDOMAINATOMS -> SIMDOMAINATOMS - External
        Total atoms in the simualtion domain
    SIMDOMAINDATA -> SIMDOMAINDATA - External
        Atomic co-ordinate  data of all the atoms
        of simualtion domain
    MAGCUBCOUNT -> COUNTER
        Counter to generate final J data
    ORGCELLREPLICA -> ORGCELLREPLICA - External
        Repetations of the original cell in all directions
    MAGCUBICNNLIST -> NNLIST - External
        Neighbor list generated by the function
    ***********************************************************/

    double XCOOD, YCOOD, ZCOOD, DISij, CX, CY, CZ, XHALF, YHALF, ZHALF, ATOMICDIS;
    XHALF = ORGCELLREPLICA*0.5;
    YHALF = XHALF;
    ZHALF = XHALF*_c2aRatio;
    //NNLIST.resize(SIMDOMAINATOMS);
    if (_c2aRatio == 1.0 ) ATOMICDIS = STRUCTCUBDIS;
    else                ATOMICDIS = STRUCTTETDIS;
    cout<<"The ATOMDIS to look for : "<<ATOMICDIS<<endl;
    for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
    {
        map <double, vector<int>> NNmap;
        for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
        {
            for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
            {
                for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                {
                    XCOOD = _SIMDOMAINDATA[i].getx();
                    YCOOD = _SIMDOMAINDATA[i].gety();
                    ZCOOD = _SIMDOMAINDATA[i].getz()*_c2aRatio;
                    if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                    else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                    if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                    else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                    if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*_c2aRatio;
                    else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*_c2aRatio;
                    for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
                    {
                        if ( j != i )
                        {
                            CX = XCOOD-_SIMDOMAINDATA[j].getx();
                            CY = YCOOD-_SIMDOMAINDATA[j].gety();
                            CZ = ZCOOD-_SIMDOMAINDATA[j].getz()*_c2aRatio;
                            DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));
                            if ( DISij <= ATOMICDIS )
                            {
                                NNmap[ATOMICDIS].push_back(j);
                            }
                        }
                    }
                }
            }
        }
        if(_c2aRatio==1.0) _SIMDOMAINDATA[i].SetSCNN(NNmap);
        else              _SIMDOMAINDATA[i].SetSTNN(NNmap);
    }
    int  NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
        {
            if ( _SIMDOMAINDATA[j].getnumber() == i )
            {

                cout<<"For Atom \""<<_atoms[i].name<<"\" Neighbours number are "<<endl;
                if(_c2aRatio==1.0){
                map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getSCNN();
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {cout<<"position of this atom: "<<_SIMDOMAINDATA[j].getx()<<" "<<_SIMDOMAINDATA[j].gety()<<" "<<_SIMDOMAINDATA[j].getz()<<endl;
                  cout << it->first << " => " << it->second.size() <<'\n';
                  vector <int> NN=it->second;
                   for (unsigned int k=0; k<NN.size();k++)
                   {
                     cout<<"number of neighbour: "<<_SIMDOMAINDATA[NN[k]].name<<" position of neighbour: "<<_SIMDOMAINDATA[NN[k]].getx()<<" "<<_SIMDOMAINDATA[NN[k]].gety()<<" "<<_SIMDOMAINDATA[NN[k]].getz()<<endl;

                   }
                 }
                 }
                 else
                 {  map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getSTNN();
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {cout<<"position of this atom: "<<_SIMDOMAINDATA[j].getx()<<" "<<_SIMDOMAINDATA[j].gety()<<" "<<_SIMDOMAINDATA[j].getz()<<endl;
                  cout << it->first << " => " << it->second.size() <<'\n';
                  vector<int> NN=it->second;
                   for (unsigned int k=0; k<NN.size();k++)
                   {
                     cout<<"number of neighbour: "<<_SIMDOMAINDATA[NN[k]].name<<" position of neighbour: "<<_SIMDOMAINDATA[NN[k]].getx()<<" "<<_SIMDOMAINDATA[NN[k]].gety()<<" "<<_SIMDOMAINDATA[NN[k]].getz()<<endl;

                   }}
                 }
                cout<<endl;
                NoCount++;
                if (NoCount > 1 ) break;
            }
        }
    }


}
void lattice::DISStrNNList(string NAME,double _c2aRatio)
{
  cout << "\nCreating "<<NAME<< " Neighbour List ... ";
    cout<<"\n";
    int ATOMi1, ATOMj1, ATOMi2, ATOMj2, DISCOL;
    double XCOOD, YCOOD, ZCOOD, DISij, ERR, CX, CY, CZ, XHALF, YHALF, ZHALF;
    XHALF = ORGCELLREPLICA*0.5;
    YHALF = XHALF;
    ZHALF = XHALF*_c2aRatio;
    if ( _c2aRatio == 1.0 )    DISCOL = 2;
    else                      DISCOL = 3;

    for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
    {   map<double,vector<int>> NNmap;
        ATOMi1 = int(_SIMDOMAINDATA[i].getnumber());

        for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
        {
            for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
            {
                for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                {
                    XCOOD = _SIMDOMAINDATA[i].getx();
                    YCOOD = _SIMDOMAINDATA[i].gety();
                    ZCOOD = _SIMDOMAINDATA[i].getz()*_c2aRatio;
                    if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                    else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                    if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                    else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                    if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*_c2aRatio;
                    else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*_c2aRatio;
                    for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
                    {
                        if ( j != i )
                        {
                            ATOMj1 = int(_SIMDOMAINDATA[j].getnumber());
                            CX = XCOOD-_SIMDOMAINDATA[j].getx();
                            CY = YCOOD-_SIMDOMAINDATA[j].gety();
                            CZ = ZCOOD-_SIMDOMAINDATA[j].getz()*_c2aRatio;
                            DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));
                            if ( DISij < MAXNNDIS )
                            {
                                for ( int VAR = 0;  VAR < STRUCTPOINTS; VAR++ )
                                {
                                    ATOMi2 = int(STRUCTDATA[VAR][0]);
                                    if ( ATOMi1 == ATOMi2 )
                                    {
                                        ATOMj2 = int(STRUCTDATA[VAR][1]);
                                        if ( ATOMj1 == ATOMj2 )
                                        {
                                            ERR = std::abs(DISij-STRUCTDATA[VAR][DISCOL]);
                                            if ( ERR < TOL )
                                            {
                                                NNmap[STRUCTDATA[VAR][DISCOL]].push_back(j);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(_c2aRatio==1.0) _SIMDOMAINDATA[i].SetSCNN(NNmap);
        else _SIMDOMAINDATA[i].SetSTNN(NNmap);
    }

    int  NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < _SIMDOMAINATOMS; j++ )
        {
            if ( _SIMDOMAINDATA[j].getnumber() == i )
            {

                cout<<"For Atom \""<<_atoms[i].name<<"\" Neighbours number are "<<endl;
                if(_c2aRatio==1.0){
                map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getSCNN();
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {cout << it->first << " => " << it->second.size() << '\n';}
                 }
                 else
                 {  map<double, vector<int>> mymap=_SIMDOMAINDATA[j].getSTNN();
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                 {cout << it->first << " => " << it->second.size() << '\n';}
                 }
                cout<<endl;
                NoCount++;
                if (NoCount > 1 ) break;
            }
        }
    }


}

int lattice::getSIMDOMAINATOMS(){
    return _SIMDOMAINATOMS;

}

atom* lattice::getSIMDOMAINDATA()
{
  return _SIMDOMAINDATA;
}

void lattice::setSIMDOMAINDATA(int i, atom thisatom)
{
  _SIMDOMAINDATA[i]=thisatom;

}


double lattice::TotalEnergy (modelstats modelnum, double kBTLnp,double gmH)
{
  //_SIGMALIST={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
/**************************************
    SIMDOMAINATOMS -> SIMDOMAINATOMS
    Total number of atoms in the sim domain
    SPIN -> SPIN
    Spin state matric of all the atoms in
    the Sim Domain
    SIGMA -> SIGMA
    Deform state matric of all the atoms in
    the Sim Domain
    MAGCUBICNNLIST -> NNLIST
    neighbor list considered for energy
    ********************************************
    Assuming only one set of structural chains..
    Needs to be modified for more than one chain
    ********************************************/

    int SPINi, SPINj, SIGMAi, SIGMAj, SIGMAj2, ATOMnn, NNChains,KSPIN;
    double HM = 0, HINT=0, HS = 0, H = 0, Jm, KnSij, SiSj, KnSSiSj;
//    KgmH=gmH*U1;
    map<double, vector<int>> mymap;
    vector <int> NNlist;
    for ( int i = 0; i < getSIMDOMAINATOMS(); i++ )
    {
        SPINi = _SPINLIST[i];
        SIGMAi = _SIGMALIST[i];
        KSPIN=_KSPINLIST[i];
        if ( SPINi != 0 )
        {
            if ( SIGMAi == 0 )                                              // Cubic State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_CubmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
                NNChains = _MCNNLIST[i].size();
                if ( NNChains > 0 )
                {    mymap=_MCNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_cubproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;}
                             else{
                              if((_subdomainnumber[i]==_subdomainnumber[ATOMnn])&(SPINi==SPINj))
                              {
                                KnSij+=1;
                              }
                             }

                          }
                          HM += -Jm*KnSij/2;

                       }
                 }

               mymap=_SCNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= 1-SIGMAj2;
                             if(_cubproblist[i]==1)
                             {if (SPINi==SPINj)
                                {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                                }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                               }
                             }
                             HINT+=Uc*2.0*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn]*KnSSiSj;
                          }
                }

                HS += - K * SiSj/2 - kBTLnp;


            }
            else                                                        // Tetragonal State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_TetmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
                //NNChains = MTNNLIST[i].size();
                NNChains = _MTNNLIST[i].size();
                if ( NNChains > 0 )
                {   mymap=_MTNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_tetproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;}
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSij+=1;
                               }
                             }

                          }
                          HM+=-Jm*KnSij/2;

                       }
                 }


             /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                mymap=_STNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= SIGMAj;
                             if(_tetproblist[i]==1)
                             {if (SPINi==SPINj)
                             {
                                 KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                             }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                  KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                               }
                             }
                             HINT+=Ut*2.0*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn]*KnSSiSj;
                          }
                }

                HS += -J * SIGMAi * SiSj/2;

                // - KgmH * SIGMAi * SiSj;
                if(_SPINLIST[i]==1) HS-=U1*gmH*SiSj*SIGMAi;
                HS-=U2*gmH*gmH*SiSj*SIGMAi/2;
            }
        }
        else
        {
            if ( SIGMAi == 0 )
            {
                SiSj = 0;
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                 mymap=_SCNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(SIGMAj==0) SiSj += 1;
                          }
                }

                HS += -K * SiSj/2 - kBTLnp;
            }
            else
            {
                SiSj = 0;

                 mymap=_STNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SIGMAj= _SIGMALIST[ATOMnn];
                            SiSj += SIGMAj;
                          }
                }
                HS += - J * SIGMAi * SiSj/2; // - KgmH * SIGMAi * SiSj;

               HS-=U2*gmH*gmH*SiSj/2*SIGMAi;
            }
        }
    }
    H =  ( HS + HM +HINT);

    //cout<<"HS "<<HS<<"HM "<<HM<<endl;
    return H;
}


double lattice::delStrEng( int sigmahext,int i, int SIGMAi, double kBTLnp,double gmH)
{

    int SIGMAj;
    double HS = 0, SIGij = 0, SIGext = 0;
    map<double, vector<int>> mymap;
    vector <int> NNlist;
    if ( SIGMAi == 0 )
    {   SIGij = 0;
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
           {NNlist=it->second;
            for(unsigned int j=0; j< NNlist.size();j++)
            {
               SIGMAj=_SIGMALIST[NNlist[j]];
               if(SIGMAj==0) SIGij+=1;
            }

           }



        HS += -2 * K * SIGij - kBTLnp;
    }
    else
    {

        SIGij = 0, SIGext = 0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
           {NNlist=it->second;
            for(unsigned int j=0; j< NNlist.size();j++)
            {
               SIGMAj=_SIGMALIST[NNlist[j]];
               SIGij+=SIGMAj;
            }

           }

        if ( SIGMAi == sigmahext )  SIGext = 1;
        HS += SIGMAi * (- 2 * J * SIGij );
        //if ( SIGMAi == 1 )  HS -= KgmH * SIGij;
        if(_SPINLIST[i]==1) HS-=U1*gmH*SIGij*SIGMAi;
        HS-=2*U2*gmH*gmH*SIGij*SIGMAi;
    }

    return HS;


}

double lattice::DeformTrialEq(int sigmahext, int i,double kBTLnp,double gmH)
{
    int SIGMAn, SIGMAi;
    double HSi=0, HSn=0, delHS=0;
    double kBT=kBTLnp/log(2);
    SIGMAi = _SIGMALIST[i];
    SIGMAn = RNDMNO.randInt(2)-1;
    if ( SIGMAi == SIGMAn ) return 0;
    HSi = delStrEng (sigmahext,i,SIGMAi,kBTLnp,gmH);
    HSn = delStrEng (sigmahext,i,SIGMAn,kBTLnp,gmH);
    delHS = HSn - HSi;
    double RNO = 0, APROB = 1;

    if ( delHS > 0 )
    {
        APROB = exp(-delHS/kBT);
        RNO = RNDMNO.randExc();
    }
    if ( RNO < APROB )
    {
        _SIGMALIST[i]=SIGMAn;
        return delHS;
    }
    else
    {
        return 0.0;
    }

}

double lattice::SpinTrial (modelstats modelnum, double kBTLnp,double gmH,int i )
{
    int SPINn, SPINi;
    double delHM,HMi,HMn;
    double kBT=kBTLnp/log(2);
    int spinhext=modelnum.getspinhext();
    SPINi = _SPINLIST[i];
    if (SPINi==0) return 0;
    int MagState = _atoms[_SIMADDRESS[i][0]].MAGSTATES;
    SPINn = RNDMNO.randInt(MagState-1)+1;
    //KSPIN[i]=MagState;
    if ( SPINi == SPINn )   return 0;

    //delHM=delTotalEnergy ( modelnum,  kBTLnp, gmH, i, SPINi,SPINn);
    delHM=delMagEng_old (spinhext, i,  SPINi,  SPINn , gmH);
    double RNO = 0, APROB = 1;
    if ( delHM > 0 )
    {
        APROB = exp(-delHM/(kBT));
        RNO = RNDMNO.randExc();
    }
    if ( RNO < APROB )
    {
        _SPINLIST[i]=SPINn;
        return delHM;
    }
    else
    {
        return 0;
    }
}

double lattice::SpinTrial_au (modelstats modelnum, double kBTLnp,double gmH,int i )
{
    int SPINn, SPINi;
    double delHM,HMi,HMn;
    double kBT=kBTLnp/log(2);
    int spinhext=modelnum.getspinhext();
    SPINi = _SPINLIST[i];
    if (SPINi==0) return 0;
    int MagState = _atoms[_SIMADDRESS[i][0]].MAGSTATES;
    SPINn = RNDMNO.randInt(MagState-1)+1;
    //KSPIN[i]=MagState;
    if ( SPINi == SPINn )   return 0;

    //delHM=delTotalEnergy_au ( modelnum,  kBTLnp, gmH, i, SPINi,SPINn);
    delHM=delMagEng_old_au (spinhext, i,  SPINi,  SPINn , gmH);
    double RNO = 0, APROB = 1;
    if ( delHM > 0 )
    {
        APROB = exp(-delHM/(kBT));
        RNO = RNDMNO.randExc();
    }
    if ( RNO < APROB )
    {
        _SPINLIST[i]=SPINn;
        return delHM;
    }
    else
    {
        return 0;
    }
}


double lattice::delMagEng_old ( int spinhext, int i, int SPINi, int SPINn ,double gmH)
{   double delHM = 0, Jm, KnSSiSj, KnSSnSj,HMi,HMn,SIGij;
    int SIGMAi, SIGMAj, NNAtoms, ATOMnn, SPINj, KnSij, KnSnj;
    SIGMAi = _SIGMALIST[i];
    map<double, vector<int>> mymap;

    if ( SIGMAi == 0 )  {
        HMi=0;HMn=0;
        if ( SPINi == spinhext )    HMi -= gmH*_CubmagmomentinSD[i];
  //      else                        HMi += gmH*_CubmagmomentinSD[i];
        if ( SPINn == spinhext )    HMn -= gmH*_CubmagmomentinSD[i];
  //      else                        HMn += gmH*_CubmagmomentinSD[i];
        if(SPINi==_KSPINLIST[i])         HMi -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
  //      else                             HMi += Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
        if(SPINn==_KSPINLIST[i])         HMn -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
  //      else                             HMn += Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
        delHM += HMn-HMi;
        //NNChains = MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0, KnSnj=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_cubproblist[i]==1)
                {if ( SIGMAj == 0 )  { if(_cubproblist[ATOMnn]==1)
                    {if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;}
                    else{
                    if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                    {if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;}
                    else{
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                    }

                    }
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{
                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                   {
                   if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                   }
                    else   {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                   }

                }

            }
            delHM -= Jm * ( KnSnj - KnSij );
        }

        /** For atoms interacting in tetra state **/
       // NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0, KnSnj = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_cubproblist[i]==1)
                {if ( SIGMAj != 0 )  {if(_tetproblist[ATOMnn]==1){
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;}
                }
                }
                else
                {
                  if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                   if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;

                  }
                }
                }

            }

            delHM -= Jm * ( KnSnj - KnSij );
        }

        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        //NNAtoms = thisatom.size();
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
 
                 if ( SIGMAj == 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 2;
                    KnSSnSj += 0;
                }
            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
             }

        //       }

       //   }
          }
        }
        delHM += Uc*2 * ( KnSSnSj - KnSSiSj )- Uc*1/2 * ( KnSnj - KnSij );

        /** Interaction Energy - Atom in NN but Tetragonal State **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
 
                 if ( SIGMAj != 0 )
                {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                }

         //      }


        //    }
          }
        }

        delHM += Ut*2 * ( KnSSnSj - KnSSiSj ) - Ut/2 * ( KnSnj - KnSij );
    }
    else    {
        HMi=0;HMn=0;
        if ( SPINi == spinhext )    HMi -= gmH*_TetmagmomentinSD[i];
 //       else                        HMi += gmH*_TetmagmomentinSD[i];
        if ( SPINn == spinhext )    HMn -= gmH*_TetmagmomentinSD[i];
 //       else                        HMn += gmH*_TetmagmomentinSD[i];
        if(SPINi==_KSPINLIST[i])         HMi -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
 //       else                             HMi += Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
        if(SPINn==_KSPINLIST[i])         HMn -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
 //       else                             HMn += Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
        delHM += HMn-HMi;
        //NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0, KnSnj=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj != 0 )  {
                          if(_tetproblist[ATOMnn]==1){
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                    }else{
                       if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                       {
                       if ( SPINi == SPINj )   KnSij += 2;
                       if ( SPINn == SPINj )   KnSnj += 2;
                       }
                       else{
                       if ( SPINi == SPINj )   KnSij += 1;
                       if ( SPINn == SPINj )   KnSnj += 1;
                       }

                    }
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{

                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                     if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                   }

                }

            }
            delHM -= Jm * ( KnSnj - KnSij );
        }


        /** For atoms interacting in cubic state **/

        //NNChains = _MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0, KnSnj = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj == 0 &&_cubproblist[ATOMnn]==1)  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{
                    if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                      if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                      }
                    }

                }


            }

            delHM -= Jm * ( KnSnj - KnSij ) ;
        }





        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0,SIGij=0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            SIGij+=SIGMAj;

                if ( SIGMAj != 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 2;
                    KnSSnSj += 0;
                }
            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }
            }





          }
        }
        HMi=0;HMn=0;
        if ( SPINi == spinhext )  HMi = U1*gmH * SIGMAi*SIGij;
        if ( SPINn == spinhext )  HMn = U1*gmH * SIGMAi*SIGij;
        delHM+=HMn-HMi;
        delHM += Ut*2 * ( KnSSnSj - KnSSiSj ) - Ut/2 * ( KnSnj - KnSij );

        /** Interaction Energy - Atom in NN but Cubic State **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];

                 if ( SIGMAj == 0 )

                 {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                 }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                 }
                }


 
          }
        }

        delHM += Uc*2 * ( KnSSnSj - KnSSiSj ) - Uc/2 * ( KnSnj - KnSij );
    }
    return delHM;
}



double lattice::delMagEng_old_au ( int spinhext, int i, int SPINi, int SPINn ,double gmH)
{   double delHM = 0, Jm, KnSSiSj, KnSSnSj,HMi,HMn,SIGij;
    int SIGMAi, SIGMAj, NNAtoms, ATOMnn, SPINj, KnSij, KnSnj;
    SIGMAi = _SIGMALIST[i];
    map<double, vector<int>> mymap;

    if ( SIGMAi == 0 )  {
        HMi=0;HMn=0;
        if ( SPINi == spinhext )    HMi -= gmH*_CubmagmomentinSD[i];
  //      else                        HMi += gmH*_CubmagmomentinSD[i];
        if ( SPINn == spinhext )    HMn -= gmH*_CubmagmomentinSD[i];
  //      else                        HMn += gmH*_CubmagmomentinSD[i];
        if(SPINi==_KSPINLIST[i])         HMi -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
  //      else                             HMi += Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
        if(SPINn==_KSPINLIST[i])         HMn -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
  //      else                             HMn += Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
        delHM += HMn-HMi;
        //NNChains = MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0, KnSnj=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_cubproblist[i]==1)
                {if ( SIGMAj == 0 )  { if(_cubproblist[ATOMnn]==1)
                    {if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;}
                    else{
                    if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                    {if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;}
                    else{
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                    }

                    }
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{
                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                   {
                   if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                   }
                    else   {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                   }

                }

            }
            delHM -= Jm * ( KnSnj - KnSij );
        }

        /** For atoms interacting in tetra state **/
       // NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0, KnSnj = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_cubproblist[i]==1)
                {if ( SIGMAj != 0 )  {if(_tetproblist[ATOMnn]==1){
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;}
                }
                }
                else
                {
                  if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                   if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;

                  }
                }
                }

            }

            delHM -= Jm * ( KnSnj - KnSij );
        }

        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        //NNAtoms = thisatom.size();
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];

                 if ( SIGMAj == 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 2;
                    KnSSnSj += 0;
                }
            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
             }

 
          }
        }
        delHM += Uc*2 * ( KnSSnSj - KnSSiSj )- Uc*1/2 * ( KnSnj - KnSij );

        /** Interaction Energy - Atom in NN but Tetragonal State **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];

                 if ( SIGMAj != 0 )
                {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }
                }

         //      }


        //    }
          }
        }

        delHM += Ut*2 * ( KnSSnSj - KnSSiSj ) - Ut/2 * ( KnSnj - KnSij );
    }
    else    {
        HMi=0;HMn=0;
        if ( SPINi == spinhext )    HMi -= gmH*_TetmagmomentinSD[i];
 //       else                        HMi += gmH*_TetmagmomentinSD[i];
        if ( SPINn == spinhext )    HMn -= gmH*_TetmagmomentinSD[i];
 //       else                        HMn += gmH*_TetmagmomentinSD[i];
        if(SPINi==_KSPINLIST[i])         HMi -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
 //       else                             HMi += Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
        if(SPINn==_KSPINLIST[i])         HMn -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
 //       else                             HMn += Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
        delHM += HMn-HMi;
        //NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0, KnSnj=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj != 0 )  {
                          if(_tetproblist[ATOMnn]==1){
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                    }else{
                       if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                       {
                       if ( SPINi == SPINj )   KnSij += 2;
                       if ( SPINn == SPINj )   KnSnj += 2;
                       }
                       else{
                       if ( SPINi == SPINj )   KnSij += 1;
                       if ( SPINn == SPINj )   KnSnj += 1;
                       }

                    }
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{

                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                     if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;
                    if ( SPINn == SPINj )   KnSnj += 2;
                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                   }

                }

            }
            delHM -= Jm * ( KnSnj - KnSij );
        }


        /** For atoms interacting in cubic state **/

        //NNChains = _MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0, KnSnj = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj == 0 &&_cubproblist[ATOMnn]==1)  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                }
                }
                else{
                    if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                      if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;
                    if ( SPINn == SPINj )   KnSnj += 1;
                      }
                    }

                }


            }

            delHM -= Jm * ( KnSnj - KnSij ) ;
        }





        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0,SIGij=0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            SIGij+=SIGMAj;

                if ( SIGMAj != 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 2;
                    KnSSnSj += 0;
                }
            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }
            }


       //       }

     //      }
          }
        }
        HMi=0;HMn=0;
        if ( SPINi == spinhext )  HMi = U1*gmH * SIGMAi*SIGij;
        if ( SPINn == spinhext )  HMn = U1*gmH * SIGMAi*SIGij;
        delHM+=HMn-HMi;
        delHM += Uc*2 * ( KnSSnSj - KnSSiSj ) - Uc/2 * ( KnSnj - KnSij );

        /** Interaction Energy - Atom in NN but Cubic State **/

        KnSij = 0, KnSSiSj = 0, KnSnj = 0, KnSSnSj = 0;
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];

                 if ( SIGMAj == 0 )

                 {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                 }
                if ( SPINn == SPINj )   {
                    //KnSnj   += 1;
                    KnSSnSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                 }
                }


       //        }

      //      }
          }
        }

        delHM += Uc*2 * ( KnSSnSj - KnSSiSj ) - Uc/2 * ( KnSnj - KnSij );
    }
    return delHM;
}



double lattice::delTotalEnergy (modelstats modelnum, double kBTLnp,double gmH,int i,int SPINi,int SPINn)
{
  //_SIGMALIST={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
/**************************************
    SIMDOMAINATOMS -> SIMDOMAINATOMS
    Total number of atoms in the sim domain
    SPIN -> SPIN
    Spin state matric of all the atoms in
    the Sim Domain
    SIGMA -> SIGMA
    Deform state matric of all the atoms in
    the Sim Domain
    MAGCUBICNNLIST -> NNLIST
    neighbor list considered for energy
    ********************************************
    Assuming only one set of structural chains..
    Needs to be modified for more than one chain
    ********************************************/

    int  SPINj, SIGMAi, SIGMAj, SIGMAj2, ATOMnn, NNChains,KSPIN;
    double HM = 0, HINT=0, HS = 0, H = 0, Jm, KnSij, SiSj, KnSSiSj;
    double HMn=0,HINTn=0,HSn=0,Hn=0,KnSnj,SnSj,KnSSnSj;
    //gmH=KgmH/U1;
    map<double, vector<int>> mymap;
    vector <int> NNlist;

        //SPINi = _SPINLIST[i];
        SIGMAi = _SIGMALIST[i];
        KSPIN=_KSPINLIST[i];
        if ( SPINi != 0 )
        {
            if ( SIGMAi == 0 )                                              // Cubic State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_CubmagmomentinSD[i];
                if ( SPINn == modelnum.getspinhext() )   HMn -= gmH*_CubmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
                if(SPINn==KSPIN)         HMn -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
                NNChains = _MCNNLIST[i].size();
                if ( NNChains > 0 )
                {    mymap=_MCNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           KnSnj=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_cubproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;
                              if(SPINn==SPINj)   KnSnj+=1;
                             }
                             else{
                              if((_subdomainnumber[i]==_subdomainnumber[ATOMnn])&(SPINi==SPINj))
                              {
                                KnSij+=1;
                              }
                              if((_subdomainnumber[i]==_subdomainnumber[ATOMnn])&(SPINn==SPINj))
                              {
                                KnSnj+=1;
                              }
                             }

                          }
                          HM += -Jm*KnSij;
                          HMn += -Jm*KnSnj;

                       }
                 }

               mymap=_SCNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0,KnSnj=0,KnSSnSj=0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             KnSSnSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= 1-SIGMAj2;
                             if(_cubproblist[i]==1)
                             {if (SPINi==SPINj)
                                {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                                }
                              if (SPINn==SPINj)
                                {
                                 KnSSnSj=0.5*(0.5-SIGMAj2)-0.25;
                                }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSSnSj=0.5*(0.5-SIGMAj2)-0.25;
                               }
                             }
                             HINT+=Uc*2.0*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn]*KnSSiSj;
                             HINTn+=Uc*2.0*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn]*KnSSnSj;
                          }
                }

                HS += - K * SiSj - kBTLnp;
                HSn += - K * SiSj - kBTLnp;


            }
            else                                                        // Tetragonal State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_TetmagmomentinSD[i];
                if ( SPINn == modelnum.getspinhext() )   HMn -= gmH*_TetmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
                if(SPINn==KSPIN)         HMn -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
                //NNChains = MTNNLIST[i].size();
                NNChains = _MTNNLIST[i].size();
                if ( NNChains > 0 )
                {   mymap=_MTNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           KnSnj=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_tetproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;
                              if(SPINn==SPINj)   KnSnj+=1;
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSij+=1;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSnj+=1;
                               }
                             }

                          }
                          HM+=-Jm*KnSij;
                          HMn+=-Jm*KnSnj;

                       }
                 }


             /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                KnSnj = 0, KnSSnSj = 0;
                mymap=_STNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= SIGMAj;
                             if(_tetproblist[i]==1)
                             {if (SPINi==SPINj)
                             {
                                 KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                             }
                             if (SPINn==SPINj)
                             {
                                 KnSSnSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                             }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                  KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                  KnSSnSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                               }
                             }
                             HINT+=Ut*2.0*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn]*KnSSiSj;
                             HINTn+=Ut*2.0*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn]*KnSSnSj;
                          }
                }

                HS += -J * SIGMAi * SiSj;
                HSn += -J * SIGMAi * SiSj;

                // - KgmH * SIGMAi * SiSj;
                if ( SPINi == modelnum.getspinhext() )  HS -= U1*gmH * SIGMAi*SiSj;
                if ( SPINn == modelnum.getspinhext() )  HSn -= U1*gmH * SIGMAi*SiSj;
                HS-=U2*gmH*gmH*SIGMAi*SiSj;
                HSn-=U2*gmH*gmH*SIGMAi*SiSj;
            }
        }
        else
        {
            if ( SIGMAi == 0 )
            {
                SiSj = 0;
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                 mymap=_SCNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(SIGMAj==0) SiSj += 1;
                          }
                }

                HS += -K * SiSj - kBTLnp;
                HSn += -K * SiSj - kBTLnp;
            }
            else
            {
                SiSj = 0;

                 mymap=_STNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SIGMAj= _SIGMALIST[ATOMnn];
                            SiSj += SIGMAj;
                          }
                }
                HS += - J * SIGMAi * SiSj; // - KgmH * SIGMAi * SiSj;
                HSn += - J * SIGMAi * SiSj;
                HS-=U2*gmH*gmH*SIGMAi*SiSj;
                HSn-=U2*gmH*gmH*SIGMAi*SiSj;
            }
        }

    H =  ( HS + HM +HINT);
    Hn=(HSn+HMn+HINTn);
    H=Hn-H;

    //cout<<"HS "<<HS<<"HM "<<HM<<endl;
    return H;
}

double lattice::delTotalEnergy_au(modelstats modelnum, double kBTLnp,double gmH,int i,int SPINi,int SPINn)
{
  //_SIGMALIST={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
/**************************************
    SIMDOMAINATOMS -> SIMDOMAINATOMS
    Total number of atoms in the sim domain
    SPIN -> SPIN
    Spin state matric of all the atoms in
    the Sim Domain
    SIGMA -> SIGMA
    Deform state matric of all the atoms in
    the Sim Domain
    MAGCUBICNNLIST -> NNLIST
    neighbor list considered for energy
    ********************************************
    Assuming only one set of structural chains..
    Needs to be modified for more than one chain
    ********************************************/

    int  SPINj, SIGMAi, SIGMAj, SIGMAj2, ATOMnn, NNChains,KSPIN;
    double HM = 0, HINT=0, HS = 0, H = 0, Jm, KnSij, SiSj, KnSSiSj;
    double HMn=0,HINTn=0,HSn=0,Hn=0,KnSnj,SnSj,KnSSnSj;
 //   gmH=KgmH/U1;
    map<double, vector<int>> mymap;
    vector <int> NNlist;

        //SPINi = _SPINLIST[i];
        SIGMAi = _SIGMALIST[i];
        KSPIN=_KSPINLIST[i];
        if ( SPINi != 0 )
        {
            if ( SIGMAi == 0 )                                              // Cubic State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_CubmagmomentinSD[i];
                if ( SPINn == modelnum.getspinhext() )   HMn -= gmH*_CubmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
                if(SPINn==KSPIN)         HMn -= -Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
                NNChains = _MCNNLIST[i].size();
                if ( NNChains > 0 )
                {    mymap=_MCNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           KnSnj=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_cubproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;
                              if(SPINn==SPINj)   KnSnj+=1;
                             }
                             else{
                              if((_subdomainnumber[i]==_subdomainnumber[ATOMnn])&(SPINi==SPINj))
                              {
                                KnSij+=1;
                              }
                              if((_subdomainnumber[i]==_subdomainnumber[ATOMnn])&(SPINn==SPINj))
                              {
                                KnSnj+=1;
                              }
                             }

                          }
                          HM += -Jm*KnSij;
                          HMn += -Jm*KnSnj;

                       }
                 }

               mymap=_SCNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0,KnSnj=0,KnSSnSj=0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             KnSSnSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= 1-SIGMAj2;
                             if(_cubproblist[i]==1)
                             {if (SPINi==SPINj)
                                {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                                }
                              if (SPINn==SPINj)
                                {
                                 KnSSnSj=0.5*(0.5-SIGMAj2)-0.25;
                                }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSSiSj=0.5*(0.5-SIGMAj2)-0.25;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSSnSj=0.5*(0.5-SIGMAj2)-0.25;
                               }
                             }
                             HINT+=Uc*2.0*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn]*KnSSiSj;
                             HINTn+=Uc*2.0*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn]*KnSSnSj;
                          }
                }

                HS += - K * SiSj - kBTLnp;
                HSn += - K * SiSj - kBTLnp;


            }
            else                                                        // Tetragonal State of SYstem
            {   if ( SPINi == modelnum.getspinhext() )   HM -= gmH*_TetmagmomentinSD[i];
                if ( SPINn == modelnum.getspinhext() )   HMn -= gmH*_TetmagmomentinSD[i];
                if(SPINi==KSPIN)         HM -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
                if(SPINn==KSPIN)         HMn -= -Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];
                //NNChains = MTNNLIST[i].size();
                NNChains = _MCNNLIST[i].size();
                if ( NNChains > 0 )
                {   mymap=_MCNNLIST[i];
                    for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                       {   Jm=it->first;
                           KnSij=0;
                           KnSnj=0;
                           NNlist=it->second;
                          for(unsigned int j=0;j<it->second.size();j++)
                          {
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(_tetproblist[i]==1)
                             {if(SPINi==SPINj)   KnSij+=1;
                              if(SPINn==SPINj)   KnSnj+=1;
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSij+=1;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                 KnSnj+=1;
                               }
                             }

                          }
                          HM+=-Jm*KnSij;
                          HMn+=-Jm*KnSnj;

                       }
                 }


             /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                KnSnj = 0, KnSSnSj = 0;
                mymap=_STNNLIST[i];
               /** Interaction Energy and lattice energy**/
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj= _SIGMALIST[ATOMnn];
                             SIGMAj2=SIGMAj*SIGMAj;
                             SiSj+= SIGMAj;
                             if(_tetproblist[i]==1)
                             {if (SPINi==SPINj)
                             {
                                 KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                             }
                             if (SPINn==SPINj)
                             {
                                 KnSSnSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                             }
                             }
                             else
                             {
                               if((SPINi==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                  KnSSiSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                               }
                               if((SPINn==SPINj)&(_subdomainnumber[i]==_subdomainnumber[ATOMnn]))
                               {
                                  KnSSnSj=-0.5 * ( 0.5 - SIGMAj2 )-0.25;
                               }
                             }
                             HINT+=Ut*2.0*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn]*KnSSiSj;
                             HINTn+=Ut*2.0*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn]*KnSSnSj;
                          }
                }

                HS += -J * SIGMAi * SiSj;
                HSn += -J * SIGMAi * SiSj;

                // - KgmH * SIGMAi * SiSj;
                if ( SPINi == modelnum.getspinhext() )  HS -= U1*gmH * SIGMAi*SiSj;
                if ( SPINn == modelnum.getspinhext() )  HSn -= U1*gmH * SIGMAi*SiSj;
                HS-=U2*gmH*gmH*SIGMAi*SiSj;
                HSn-=U2*gmH*gmH*SIGMAi*SiSj;
            }
        }
        else
        {
            if ( SIGMAi == 0 )
            {
                SiSj = 0;
                SiSj = 0, KnSij = 0, KnSSiSj = 0;
                 mymap=_SCNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SPINj=_SPINLIST[ATOMnn];
                             SIGMAj=_SIGMALIST[ATOMnn];
                             if(SIGMAj==0) SiSj += 1;
                          }
                }

                HS += -K * SiSj - kBTLnp;
                HSn += -K * SiSj - kBTLnp;
            }
            else
            {
                SiSj = 0;

                 mymap=_STNNLIST[i];
                for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
                {
                    NNlist=it->second;
                    for(unsigned int j=0;j<it->second.size();j++)
                          {  KnSSiSj=0;
                             ATOMnn=NNlist[j];
                             SIGMAj= _SIGMALIST[ATOMnn];
                            SiSj += SIGMAj;
                          }
                }
                HS += - J * SIGMAi * SiSj; // - KgmH * SIGMAi * SiSj;
                HSn += - J * SIGMAi * SiSj;
                HS-=U2*gmH*gmH*SIGMAi*SiSj;
                HSn-=U2*gmH*gmH*SIGMAi*SiSj;
            }
        }

    H =  ( HS + HM +HINT);
    Hn=(HSn+HMn+HINTn);
    H=Hn-H;

    //cout<<"HS "<<HS<<"HM "<<HM<<endl;
    return H;
}



double lattice::DeformTrialPd (int sigmahext,int spinhext,int i,double kBT,double kBTLnp,double gmH)
{
    int SIGMAn, SIGMAi,SPINi;
    double HSi, HSn, delHS;

    SIGMAi = _SIGMALIST[i];
    SIGMAn = RNDMNO.randInt(2)-1;
    if ( SIGMAi == SIGMAn ) return 0;
    HSi = delStrEng (sigmahext,i,SIGMAi,kBTLnp,gmH);
    HSn = delStrEng (sigmahext,i,SIGMAn,kBTLnp,gmH);
    delHS = HSn - HSi;
    double RNO = 0, APROB = 1;
    if ( delHS > 0 )
    {
        APROB = exp(-delHS/kBT);
        RNO = RNDMNO.randExc();
    }
    if ( RNO < APROB )
    {
        SPINi = _SPINLIST[i];
        if ( SPINi != 0 ) delHS += ( IntMagEng(spinhext,i, SPINi, SIGMAn,gmH) - IntMagEng(spinhext,i, SPINi, SIGMAi,gmH) );
        _SIGMALIST[i]=SIGMAn;
        return delHS;
    }
    else
    {
        return 0.0;
    }
}

double lattice::IntMagEng ( int spinhext, int i, int SPINi, int SIGMAi ,double gmH)
{
    double delHM = 0, Jm, KnSSiSj, HMi;
    int SIGMAj, NNAtoms, ATOMnn, SPINj, KnSij;
    map<double, vector<int>> mymap;

    if ( SIGMAi == 0 )  {
        HMi=0;//HMn=0;
        if ( SPINi == spinhext )    HMi -= gmH*_CubmagmomentinSD[i];
        else                        HMi -= 0;
        if(SPINi==_KSPINLIST[i])         HMi -= Kaniscub*_CubmagmomentinSD[i]*_CubmagmomentinSD[i];
//        NNChains = _MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_cubproblist[i]==1){
                if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;

                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;

                }
                }
                else{
                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
                   {
                     if ( SIGMAj == 0 )  {
                     if ( SPINi == SPINj )   KnSij += 2;

                      }
                    else    {
                       if ( SPINi == SPINj )   KnSij += 1;

                      }


                   }

                }

            }
            delHM -= Jm * ( KnSij ) ;
        }


        /** For atoms interacting in tetra state **/
//        NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;

                }
                }
                else{
                    if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                      if ( SIGMAj != 0 )  {
                      if ( SPINi == SPINj )   KnSij += 1;

                      }

                    }


                }


            }

            delHM -= Jm * (  KnSij ) ;
        }

        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0;
        //NNAtoms = thisatom.size();
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            if(_cubproblist[i]==1){
            if ( SIGMAj == 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }

            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }

            }
            }
            else{
               if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
             if ( SIGMAj == 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }

            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }

            }


               }


            }
          }
        }
        delHM += Uc*2.0 * (  KnSSiSj ) - Uc* 0.5 * ( KnSij );

        /** Interaction Energy - Atom in NN but Tetragonal State **/

        KnSij = 0, KnSSiSj = 0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            if(_cubproblist[i]==1){
            if ( SIGMAj != 0 )

            {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }

            }
            }
            else{

               if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
            if ( SIGMAj != 0 )

            {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_CubmagmomentinSD[i]*_TetmagmomentinSD[ATOMnn];
                }

            }


               }


            }
          }
        }

        delHM += Ut*2.0 * (  KnSSiSj ) - Ut*0.5 * (  KnSij );
    }
    else    {
        HMi=0;
        if ( SPINi == spinhext )    HMi -= gmH*_TetmagmomentinSD[i];
        else                        HMi -= 0;
        if(SPINi==_KSPINLIST[i])         HMi -= Kanistet*_TetmagmomentinSD[i]*_TetmagmomentinSD[i];

//        NNChains = _MTNNLIST[i].size();
        mymap=_MTNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
        {   Jm=it->first;
            NNAtoms=it->second.size();
            KnSij=0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;

                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;

                }
                }
                else{
                   if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                   if ( SIGMAj != 0 )  {
                    if ( SPINi == SPINj )   KnSij += 2;

                }
                else    {
                    if ( SPINi == SPINj )   KnSij += 1;

                }


                   }

                }

            }
            delHM -= Jm * (  KnSij ) ;
        }


        /** For atoms interacting in cubic state **/

//        NNChains = _MCNNLIST[i].size();
        mymap=_MCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
            Jm = it->first;
            NNAtoms=it->second.size();
            KnSij = 0;
            for (int k=0;k<NNAtoms;k++){
                ATOMnn = it->second[k];
                SPINj = _SPINLIST[ATOMnn];
                SIGMAj = _SIGMALIST[ATOMnn];
                if(_tetproblist[i]==1){
                if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;

                }
                }
                else{
                  if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
                    if ( SIGMAj == 0 )  {
                    if ( SPINi == SPINj )   KnSij += 1;

                }


                  }

                }


            }

            delHM -= Jm * ( KnSij ) ;
        }





        /** Interaction Energy **/

        KnSij = 0, KnSSiSj = 0;
        mymap=_STNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            if(_tetproblist[i]==1){
            if ( SIGMAj != 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }

            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }

            }
            }
            else{

              if(_subdomainnumber[i]==_subdomainnumber[ATOMnn]){
             if ( SIGMAj != 0 )  {
                if ( SPINi == SPINj )   {
                    //KnSij   += 2;
                    KnSSiSj += 0;
                }

            }
            else    {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }

            }

              }

            }
          }
        }

        delHM += Ut*2.0 * (  KnSSiSj ) - Ut*0.5 * (  KnSij );

        /** Interaction Energy - Atom in NN but Cubic State **/

        KnSij = 0, KnSSiSj = 0;
        mymap=_SCNNLIST[i];
        for (std::map<double,vector<int>>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
        NNAtoms=it->second.size();
        for ( int j = 0; j < NNAtoms; j++ ) {

            ATOMnn = it->second[j];
            SIGMAj = _SIGMALIST[ATOMnn];
            SPINj  = _SPINLIST[ATOMnn];
            if(_tetproblist[i]==1){
            if ( SIGMAj == 0 )

            {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }

            }
            }
            else{
              if(_subdomainnumber[i]==_subdomainnumber[ATOMnn])
              {
            if ( SIGMAj == 0 )

            {
                if ( SPINi == SPINj )   {
                    //KnSij   += 1;
                    KnSSiSj -= 0.5*_TetmagmomentinSD[i]*_CubmagmomentinSD[ATOMnn];
                }

            }

              }

            }
          }
        }

        delHM += Uc*2.0 * (  KnSSiSj ) - Uc*0.5 * (  KnSij );
    }
    return delHM;
}

double lattice::NMCalcnew ()
{  // _SPINLIST={0,0,0,0,1,2,2,4,1,1,1,2,2,2,3,3};
    double MSNO, mcode;
    double NM = 0;

    vector<double>NCOUNT;
    vector<double>NMAX;
    mcode=0;
    for(int i=0;i<EAATOMTOMS;i++){
     if(mcode<_atoms[i].MAGCODES) mcode=_atoms[i].MAGCODES;

    }
    std::vector<vector<int> > NSTATES(mcode);
    NCOUNT.resize(mcode);
    NMAX.resize(mcode);
    for(int i=0;i<mcode;i++) NCOUNT[i]=0;
    for(int i=0;i<EAATOMTOMS;i++){
      NSTATES[_atoms[i].MAGCODES-1].resize(_atoms[i].MAGSTATES);
      NCOUNT[_atoms[i].MAGCODES-1]+=_atoms[i].getcount();
    }
     //cout<<NSTATES[1][2];
    for (int i=0; i<_SIMDOMAINATOMS;i++)
    {  MSNO=_SIMDOMAINDATA[i].MAGSTATES;
       if(MSNO>0)
       {
         NSTATES[_SIMDOMAINDATA[i].MAGCODES-1][_SPINLIST[i]-1]+=1;

       }

    }
    for(int i=0;i<mcode;i++)
    {
      for(unsigned int j=0;j<NSTATES[i].size();j++){
         if(NMAX[i]<NSTATES[i][j]) NMAX[i]=NSTATES[i][j];
      }
      if(NMAX[i]>0){
      NM+=(NSTATES[i].size()*NMAX[i]-NCOUNT[i])/(NSTATES[i].size()-1);

      }
    }


    return NM;
}

double lattice::NMCalc ()
{  //_SPINLIST={0,0,0,0,1,2,3,4,1,1,1,1,2,2,3,3};
   int MSNO, BEFORE, j, NCOUNT, MAGXT, NMAX, MAGC;
    double NM = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        MSNO = _atoms[i].MAGSTATES;
        MAGC = _atoms[i].MAGCODES;
        if ( MSNO > 0 )
        {
            j = i - 1;
            BEFORE = 0;
            while ( j >= 0 )
            {
                if ( MAGC == _atoms[j].MAGCODES )
                {
                    BEFORE = 1;
                    break;
                }
                j--;
            }
            if ( ! BEFORE )
            {
                int *NSTATES;
                NSTATES = new int [MSNO];
                NCOUNT = 0;
                for ( int k = 0; k < MSNO; k++ )
                {
                    NSTATES[k] = 0;
                }
                for ( int k = 0; k < _SIMDOMAINATOMS; k++ )
                {
                    if ( int(_SIMADDRESS[k][0]) == i )
                    {
                        NSTATES[_SPINLIST[k]-1] += 1;
                    }
                }
                NCOUNT += _atoms[i].getcount();
                for ( int j = (i+1); j < EAATOMTOMS; j++ )
                {
                    MAGXT = _atoms[j].MAGCODES;
                    if ( MAGXT == MAGC )
                    {
                        for ( int k = 0; k < _SIMDOMAINATOMS; k++ )
                        {
                            if ( int(_SIMADDRESS[k][0]) == j )
                            {
                                NSTATES[_SPINLIST[k]-1] += 1;
                            }
                        }
                        NCOUNT += _atoms[j].getcount();
                    }
                }
                NMAX = NSTATES[0];
                for (int MAX = 1; MAX < MSNO; MAX++)
                {
                    if ( NSTATES[MAX] > NMAX)
                    {
                        NMAX = NSTATES[MAX];
                    }
                }
                NM += ( MSNO*NMAX - NCOUNT ) / ( MSNO - 1 );
                delete [] NSTATES;
            }
        }
    }
    return NM;
}
double lattice::TDCalc ()
{
    double EPS = 0;
    for ( int i = 0; i < _SIMDOMAINATOMS; i++ )
    {
        EPS +=_SIGMALIST[i];
    }
    EPS = EPS/_SIMDOMAINATOMS;
    return EPS;
}

void lattice::convertSIMDOMAINDATA()
{
   for(int i=0;i<_SIMDOMAINATOMS;i++)
   {
     vector< double > oneatom;
     oneatom.push_back(_SIMDOMAINDATA[i].getnumber());
     oneatom.push_back(_SIMDOMAINDATA[i].getx());
     oneatom.push_back(_SIMDOMAINDATA[i].gety());
     oneatom.push_back(_SIMDOMAINDATA[i].getz());
     _SIMADDRESS.push_back(oneatom);
     int spin=_SIMDOMAINDATA[i].getspin();
     _SPINLIST.push_back(spin);
     int sigma=_SIMDOMAINDATA[i].getsigma();
     _SIGMALIST.push_back(sigma);
     int kspin=_SIMDOMAINDATA[i].getkspin();
     _KSPINLIST.push_back(kspin);
     _MCNNLIST.push_back(_SIMDOMAINDATA[i].getMCNN());
     _MTNNLIST.push_back(_SIMDOMAINDATA[i].getMTNN());
     _SCNNLIST.push_back(_SIMDOMAINDATA[i].getSCNN());
     _STNNLIST.push_back(_SIMDOMAINDATA[i].getSTNN());
     _CubmagmomentinSD.push_back(_SIMDOMAINDATA[i].CubMagMoment);
     _TetmagmomentinSD.push_back(_SIMDOMAINDATA[i].TetMagMoment);
   }


}


void lattice::creatsubdomain()
{   double l=xsub;
    double j=ysub;
    double k=zsub;
    double x1,y1,z1;//dimension of one subdomain
    x1=double(ORGCELLREPLICA/l);
    y1=double(ORGCELLREPLICA/j);
    z1=double(ORGCELLREPLICA/k);
    vector <vector < int > > subdomainlist(l*j*k);
    int r=0;
    for(int i=0;i<l;i++){
       for(int ii=0;ii<j;ii++){
          for(int iii=0;iii<k;iii++){
              subdomainlist[r]={i,ii,iii};
              r++;

          }
       }

    }
    map<vector<int>,int> subdomainmap;
    for (unsigned int i=0;i<subdomainlist.size();i++){
      subdomainmap[subdomainlist[i]]=i;

    }
    //map<int,vector<int > > _subdomainnnlist;
   for(unsigned int i=0;i<subdomainlist.size();i++){

       for (unsigned int ii=0;ii<subdomainlist.size();ii++){
           if(i!=ii)
           {  int x,y,z,x2,y2,z2;
               x=subdomainlist[i][0];
               y=subdomainlist[i][1];
               z=subdomainlist[i][2];
               x2=subdomainlist[ii][0];
               y2=subdomainlist[ii][1];
               z2=subdomainlist[ii][2];
              if((abs(x-x2)==1&&y==y2&&z==z2)||(abs(y-y2)==1&&x==x2&&z==z2)||(abs(z-z2)==1&&x==x2&&y==y2))
               {
                 _subdomainnnlist[i].push_back(ii);
               }
              if(subdomainlist[i][0]==0||subdomainlist[i][0]==xsub-1)
              { if(subdomainlist[i][1]==subdomainlist[ii][1]&&subdomainlist[i][2]==subdomainlist[ii][2])
                  {
                   if(subdomainlist[ii][0]==l-1-x) _subdomainnnlist[i].push_back(ii);
                  }

              }
              if(subdomainlist[i][1]==0||subdomainlist[i][1]==ysub-1)
              { if(subdomainlist[i][0]==subdomainlist[ii][0]&&subdomainlist[i][2]==subdomainlist[ii][2])
                  {
                    if(subdomainlist[ii][1]==j-1-y) _subdomainnnlist[i].push_back(ii);
                  }

              }
              if(subdomainlist[i][2]==0||subdomainlist[i][2]==zsub-1)
              { if(subdomainlist[i][1]==subdomainlist[ii][1]&&subdomainlist[i][0]==subdomainlist[ii][0])
                  {
                    if(subdomainlist[ii][2]==k-1-z) _subdomainnnlist[i].push_back(ii);
                  }

              }
           }


       }

   }
   //vector<int> _subdomainnumber;
    for(int i=0;i<_SIMDOMAINATOMS;i++){
        int xx,yy,zz;
        xx=int(_SIMADDRESS[i][1]/x1);
        yy=int(_SIMADDRESS[i][2]/y1);
        zz=int(_SIMADDRESS[i][3]/z1);
        _subdomainnumber.push_back(subdomainmap[{xx,yy,zz}]);
    }


}

void lattice::creatprob(double gmH)
{
   double prob1,prob2,RNO=0;
   for(int i=0;i<_SIMDOMAINATOMS;i++)

     { atom at1=_SIMDOMAINDATA[i];
       prob1=exp(-(SQUARE(_CubmagmomentinSD[i])*std::abs(Kaniscub))/(std::abs(gmH)*std::abs(_CubmagmomentinSD[i])));
       prob2=exp(-((_TetmagmomentinSD[i]*_TetmagmomentinSD[i])*std::abs(Kanistet))/(std::abs(gmH)*std::abs(_TetmagmomentinSD[i])));
       RNO=RNDMNO.rand();
       if(POLY!=0){
       if(RNO<=prob1){
          _cubproblist.push_back(1);
       }
       else{
         _cubproblist.push_back(0);
       }
       if(RNO<=prob2){
          _tetproblist.push_back(1);
       }
       else{
          _tetproblist.push_back(0);
       }
       }
       else{
          _cubproblist.push_back(1);
          _tetproblist.push_back(1);

       }
     }


}

void lattice::assignkspin()
{
   int n=xsub*ysub*zsub;
   int m=100;
   vector <map<int,int> > mlist;


   for(int i=0;i<n;i++)

   { m=RNDMNO.randInt(3-1)+1;
     map <int,int> mymap;
    for(int j=0;j<EAATOMTOMS;j++)
   {
      if(_atoms[j].MAGSTATES>0) mymap[j]=m;
      else mymap[j]=0;

   }

     mlist.push_back(mymap);

   }
   for(int i=0; i<_SIMDOMAINATOMS;i++){
     _SIMDOMAINDATA[i].setkspin(mlist[_subdomainnumber[i]][_SIMDOMAINDATA[i].getnumber()]);
     _KSPINLIST[i]=mlist[_subdomainnumber[i]][_SIMDOMAINDATA[i].getnumber()];

   }



}

double lattice::MagCalc( int isau)
{ int MSNO;
  double  MA=0;
  double MM=0;
  double cubfrac=0;
  double tetfrac=0;

    vector<double>NCOUNT;
    vector<double>NMAX;
    std::vector<vector<int> > NSTATES(EAATOMTOMS);
    NCOUNT.resize(EAATOMTOMS);
    NMAX.resize(EAATOMTOMS);
    for(int i=0;i<EAATOMTOMS;i++) {
    NCOUNT[i]=0;
    NMAX[i]=0;
    }
    for(int i=0;i<EAATOMTOMS;i++){
      NSTATES[i].resize(_atoms[i].MAGSTATES);
      NCOUNT[i]=_atoms[i].getcount();
    }
    for (int i=0; i<_SIMDOMAINATOMS;i++)
    {  MSNO=_SIMDOMAINDATA[i].MAGSTATES;
       if(MSNO>0)
       {
         NSTATES[_SIMDOMAINDATA[i].getnumber()][_SPINLIST[i]-1]+=1;

       }


    }
    for(int i=0;i<EAATOMTOMS;i++)
    {
      for(unsigned int j=0;j<NSTATES[i].size();j++){
         if(NMAX[i]<NSTATES[i][j]) NMAX[i]=NSTATES[i][j];
      }

      _atoms[i].setmag((NSTATES[i].size()*NMAX[i]-NCOUNT[i])/(NSTATES[i].size()-1)/NCOUNT[i]);


    }


 
    for(int i=0;i<_SIMDOMAINATOMS;i++)
    {
       if(_SIGMALIST[i]==0) cubfrac+=1;
       else tetfrac+=1;

    }

    cubfrac=cubfrac/_SIMDOMAINATOMS;
    tetfrac=tetfrac/_SIMDOMAINATOMS;
    for(int i=0;i<EAATOMTOMS;i++){
    MA+=_atoms[i].CubMagMoment*_atoms[i].getmag()*NCOUNT[i]/(_SIMDOMAINATOMS/4);
    MM+=_atoms[i].TetMagMoment*_atoms[i].getmag()*NCOUNT[i]/(_SIMDOMAINATOMS/4);
    }
    return cubfrac*MA+tetfrac*MM;
}

void lattice::assignspin()
{
   for(int i=0;i<_SIMDOMAINATOMS;i++){
      if(_SPINLIST[i]!=0)   _SPINLIST[i]=_KSPINLIST[i];
   }

}

void lattice::assignsigma()
{
   for(int i=0; i<_SIMDOMAINATOMS; i++){
      _SIGMALIST[i]=RNDMNO.randInt(2)-1;

   }

}
