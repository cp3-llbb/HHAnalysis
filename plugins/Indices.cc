#include <string>

#include <cp3_llbb/HHAnalysis/interface/Indices.h>

namespace HHAnalysis {
  
  // Combination of lepton ID + lepton Isolation for a single lepton
  uint16_t lepIDIso(const lepID::lepID& id, const lepIso::lepIso& iso){
    return lepIso::Count * id + iso;
  }
  std::string lepIDIsoStr(const lepID::lepID& id, const lepIso::lepIso& iso){
    return "ID" + lepID::map.at(id) + "_Iso" + lepIso::map.at(iso);
  }

  // Combination of lepton ID + lepton Isolation for a Dilepton object
  uint16_t leplepIDIso(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2){
    return lepIso::Count*lepID::Count*lepIso::Count * id1 + lepID::Count*lepIso::Count * iso1 + lepIso::Count * id2 + iso2;
  }
  std::string leplepIDIsoStr(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2){
    return "ID" + lepID::map.at(id1) + lepID::map.at(id2) + "_Iso" + lepIso::map.at(iso1) + lepIso::map.at(iso2);
  }

  // Combination of jet ID and B-tagging working point for one jet (not yet available)
  uint16_t jetIDbtagWP(const jetID::jetID& id, const btagWP::btagWP& wp){
    return btagWP::Count * id + wp;
  }
  std::string jetIDbtagWPStr(const jetID::jetID& id, const btagWP::btagWP& wp){
    return "ID" + jetID::map.at(id) + "_B" + btagWP::map.at(wp);
  }
  
  // Combination of jet ID and B-tagging working points for two jets 
  uint16_t jetjetIDbtagWP(const jetID::jetID& id1, const btagWP::btagWP& wp1, const jetID::jetID& id2, const btagWP::btagWP& wp2){
    return 
      btagWP::Count*jetID::Count*btagWP::Count * id1 + 
                 jetID::Count*btagWP::Count * wp1 + 
                              btagWP::Count * id2 + 
                                           wp2 ;
  }
  std::string jetjetIDbtagWPStr(const jetID::jetID& id1, const jetID::jetID& id2, const btagWP::btagWP wp1, const btagWP::btagWP wp2){
    return "ID" + jetID::map.at(id1) + jetID::map.at(id2) + "_B" + btagWP::map.at(wp1) + btagWP::map.at(wp2);
  }
  
  // Combination of lepton ID, lepton Isolation, jet ID, B-tagging working points and jetPair ordering for a two-lepton-two-b-jets object
  uint16_t leplepIDIsojetjetIDbtagWPPair(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2, const jetID::jetID& jetid1, const jetID::jetID& jetid2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2, const jetPair::jetPair& jetpair){
    return 
      lepIso::Count*lepID::Count*lepIso::Count*jetID::Count*jetID::Count*btagWP::Count*btagWP::Count*jetPair::Count * id1     + 
                    lepID::Count*lepIso::Count*jetID::Count*jetID::Count*btagWP::Count*btagWP::Count*jetPair::Count * iso1    + 
                                 lepIso::Count*jetID::Count*jetID::Count*btagWP::Count*btagWP::Count*jetPair::Count * id2     + 
                                               jetID::Count*jetID::Count*btagWP::Count*btagWP::Count*jetPair::Count * iso2    + 
                                                            jetID::Count*btagWP::Count*btagWP::Count*jetPair::Count * jetid1  + 
                                                                         btagWP::Count*btagWP::Count*jetPair::Count * jetid2  + 
                                                                                       btagWP::Count*jetPair::Count * wp1     + 
                                                                                                     jetPair::Count * wp2     + 
                                                                                                                      jetpair ; 
  }
  std::string leplepIDIsojetjetIDbtagWPPairStr(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2, const jetID::jetID& jetid1, const jetID::jetID& jetid2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2, const jetPair::jetPair& jetpair){
    return "lep_ID" + lepID::map.at(id1) + lepID::map.at(id2) + "_Iso" + lepIso::map.at(iso1) + lepIso::map.at(iso2) + "_jet_ID" + jetID::map.at(jetid2) + jetID::map.at(jetid2) + "_B" + btagWP::map.at(wp1) + btagWP::map.at(wp2) + "_Ordered" + jetPair::map.at(jetpair);
  }

}
