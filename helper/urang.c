double uran(long* iudum) {
  int IA=16807;
  long IM=2147483647;
  double AM=(1.0/IM);
  int IQ=127773;
  int IR=2836;
  int NTAB=32;
  double NDIV=(1+(IM-1)/NTAB);
  double RNMX=(1.0-1.2e-7);
  int j;
  long k;
  static long iy=0;
  static long iv[32];
  double temp;

  if(*iudum <=0 || !iy) {
    if(-(*iudum) < 1) *iudum=1;
    else *iudum = -(*iudum);
    for(j=NTAB+7;j>=0;j--) {
      k=*iudum/IQ;
      *iudum=IA*(*iudum-k*IQ)-IR*k;
      if(*iudum < 0) *iudum = *iudum + IM;
      if(j<NTAB) iv[j] = *iudum;
    }
    iy=iv[0];
  }
  k=*iudum/IQ;
  *iudum=IA*(*iudum-k*IQ)-IR*k;
  if(*iudum<0) *iudum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *iudum;
  if((temp=AM*iy) > RNMX) return RNMX;
  else {
   return temp;
  }
}
