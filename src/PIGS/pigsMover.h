#ifndef PIGSMOVER_H
#define PIGSMOVER_H

#include "PIGSpropagator.h"
class towerSampling;

class pigsMove
{
public:
  typedef qmcDriver_t::configurations_t configurations_t;
  typedef system_t::rand_t rand_t;
  typedef pairParticleApproximationChain propagator_t;
  
  pigsMove(){nAttempts=0;nSuccess=0;}
  
  virtual void move(configurations_t & configurations,propagator_t &prop)=0;
  
  double getAcceptanceRatio(){return nSuccess/nAttempts;}
  
  rand_t * getRandomGenerator(){return randO;}
  rand_t * setRandomGenerator(rand_t * randO2){randO=randO2;}
  
  virtual string getName(){return "pigsMove";}
protected:
  void recordSuccessfullMove(){nAttempts+=1;nSuccess+=1;}
  void recordUnSuccessfullMove(){nAttempts+=1;}
  
private:
  double nAttempts;
  double nSuccess;
  rand_t * randO;
};

class pigsMover
  {
    typedef typename system_t::rand_t rand_t;
    typedef qmcDriver_t::propagator_t propagator_t;
    typedef typename qmcDriver_t::configurations_t configurations_t;
    typedef typename system_t::pos_t pos_t;
    typedef pigsMove move_t;
    
  public:
    
    pigsMover();
    void move(configurations_t & configurations);
    void unmove(configurations_t & configurations);
    void addMove(move_t & move,double weight);
    void setPropagator(propagator_t * prop){propagator=prop;};
    propagator_t * getPropagator(){return propagator;}
    void setRandomGenerator(rand_t* rand2){randO=rand2;};
    string info();
    void setSampler(towerSampling* sampler2){sampler=sampler2;} ;
    
  private:
    rand_t *randO;
    propagator_t *propagator;
    vector<move_t *> movesToAttempt;
    vector<double> movesWeights;
    towerSampling *sampler;
};

class levyReconstructor
{
  typedef system_t::rand_t rand_t;
  typedef qmcDriver_t::configurations_t configurations_t;
  typedef system_t::geometry_t geometry_t;
  typedef system_t::pos_t pos_t;
public:
  
  void save(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,const configurations_t & configurations);
  void construct(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t &configurations);
  void recover(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations);
  
  void setRandomGenerator(rand_t * randO2){randO=randO2;}
  void setGeometry(geometry_t* geo_){geo=geo_;}
  geometry_t* getGeometry(){return geo;}
  void setConfigurationsBackup(configurations_t * configurations2){configurationsBackup=configurations2;gaussianVariables.resize(configurations2->size());};
  
  void setTimeStep(double timeStep_){timeStep=timeStep_;}
  
 
private:
  
  geometry_t* geo;
  rand_t* randO;
  double timeStep;
  vector<double> gaussianVariables;
  configurations_t * configurationsBackup;
};

class wiggleMove : public pigsMove
{
  
public:
  virtual string getName(){return "wiggle";}
  wiggleMove(){iParticle=0;i1=0;i2=0;};
  virtual void move(configurations_t & configurations,propagator_t & prop);
  
  void setConstructor(levyReconstructor* reconstructor_){reconstructor=reconstructor_;}
  
private:
  levyReconstructor *reconstructor;
  int iParticle;int i1;int i2;
};


class moveVariationalTails : public pigsMove
{
public:
  typedef system_t::pos_t pos_t;
  typedef system_t::geometry_t geometry_t;
  string getName(){return "tails";}
  moveVariationalTails(){gaussianRandomNumbers.resize(4);};
  virtual void move(configurations_t & configurations,propagator_t & prop);
  void setGeometry(geometry_t * geo_){geo=geo_;}
  void setTimeStep(double timeStep_){timeStep=timeStep_;}
private:
  vector<double> gaussianRandomNumbers;
  double timeStep;
  geometry_t * geo;
};


class translateMove : public pigsMove

{
public:
  
  translateMove(){displacement=0.01;}
  virtual string getName(){return "trans.";}
  void setDisplacement(double dis){displacement=0.1;}
  void setConstructor(levyReconstructor* reconstructor_){reconstructor=reconstructor_;}
  
  virtual void move(configurations_t & configurations,propagator_t & prop);
  
private:
  levyReconstructor *reconstructor;
  double displacement;
  
};

class moveSingleBead : public pigsMove
{
  typedef system_t::geometry_t geometry_t;
public:
  moveSingleBead(){displacement=0.01;};
  virtual string getName(){return "singleBead";}
  void setDisplacement(double dis){displacement=dis;}
  
  virtual void move(configurations_t & configurations,propagator_t & prop);
  void setGeometry(geometry_t * geo2){geo=geo2;}
  
private:
  geometry_t * geo;
  double displacement;
};
  
#endif
