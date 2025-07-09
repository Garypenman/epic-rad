#pragma once
/**
 * Derived class of ParticleModifier to add ePIC specific utilities
 */

#include "ParticleModifier.h"
#include <TMath.h>

namespace rad{
  namespace epic{


    ///////////////////////////////////////////////////////////////
    ///Class Definition
    ///////////////////////////////////////////////////////////////
    class ePICParticleModifier : public rad::config::ParticleModifier {
      
    public:
      
    ePICParticleModifier(rad::config::ConfigReaction& cr): ParticleModifier{cr}{};
 
     


    };




    
  }//epic
}//rad
