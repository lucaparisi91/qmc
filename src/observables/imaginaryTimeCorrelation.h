#ifndef IMAGINARYTIMECORRELATION_H
#define IMAGINARYTIMECORRELATION_H


void timeDifferenceSquaresAverage(const vector<double> & positions,vector<double> & timeCorrelations);

void timeProductAverage(const vector<double> & positions,vector<double> & timeCorrelations);


template<class walker_t,class wave_t>
class windingNumber : public measurement<walker_t,wave_t,measure_vector >
{
public:
  windingNumber(measure_vector* ms_,int bins_,int setA_,int indexStorage_) : measurement<walker_t,wave_t,measure_vector >(ms_){indexStorage=indexStorage_;bins=bins_;setA=setA_;timeCorrelationsCurrent.resize(bins);};
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    // store the measurements on the walker
    
    w->md[indexStorage]->add(centerOfMassNoBC((*w->state)[setA]));
    
    // if the time of the windows are espired
    
    if (w->md[indexStorage]->isFilled()==1)
      {
	timeDifferenceSquaresAverage( w->md[indexStorage]->currentVector(),timeCorrelationsCurrent);
	w->md[indexStorage]->reset();
	this->ms->add(timeCorrelationsCurrent,0);
      }
    
  }
  
private:
  int bins;
  vector<double> timeCorrelationsCurrent;
  int indexStorage;
  int setA;
};

template<class walker_t,class wave_t>
class windingNumberSpin : public measurement<walker_t,wave_t,measure_vector >
{
public:
  windingNumberSpin(measure_vector* ms_,int bins_,int setA_,int indexStorage_) : measurement<walker_t,wave_t,measure_vector >(ms_){indexStorage=indexStorage_;bins=bins_;setA=setA_;timeCorrelationsCurrent.resize(bins);};
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    // store the measurements on the walker
    
    w->md[indexStorage]->add(centerOfMassSpinNoBC((*w->state)[setA]));
    
    // if the time of the windows are espired
    
    if (w->md[indexStorage]->isFilled()==1)
      {
	timeDifferenceSquaresAverage( w->md[indexStorage]->currentVector(),timeCorrelationsCurrent);
	w->md[indexStorage]->reset();
	this->ms->add(timeCorrelationsCurrent,0);
      }
    
  }
  
private:
  int bins;
  vector<double> timeCorrelationsCurrent;
  int indexStorage;
  int setA;
};

#endif
