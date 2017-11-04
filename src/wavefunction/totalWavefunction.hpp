template< class tm>
void total_wavefunction<tm>::laplacianGradient(const all_particles_t & p,value_t & e,value_t & eF, grad_t & grad)
{
  value_t tmp;
  int i,j;
  
  grad.reset();
  
  e=0;
  tmp=0;
  for(i=0;i<waves.size();i++)
    {
      waves[i]->laplacianMinusGradientSquared(p,grad,tmp);
      e=e+tmp;
    }
  //cout << "K:" <<e<<endl;
  eF=0;
  
  // compute the drift force term
  for(i=0;i<grad.size();i++)
    {
      for(j=0;j<grad[i].size();j++)
	{
	  eF+=pow(grad[i][j],2);
	}
    }
  //cout<<"V: "<<e_f<<endl;
  // sum the two contributions
  e=e + eF;
  // multiply for the diffusion coefficient
  e=-e/2.;
  eF=eF/2;
}

template< class tm>
void total_wavefunction<tm>::gradient(const all_particles_t & p,grad_t & grad)
{
  value_t tmp;
  int i,j;
  
  grad.reset();
  
  for(i=0;i<waves.size();i++)
    {
      waves[i]->gradient(p,grad);
    }
}
