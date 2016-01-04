#pragma once

#include <map>
#include <array>
#include <string>

namespace HHAnalysis {
  
  // lepton IDs
  namespace lepID {
    enum lepID{ L, T, Count };
    // Ugly way to allow iterating over all items in the enumeration ( for(const lepID::lepID& id: lepID::it) )
    const std::array<lepID, Count> it = {{ L, T }};
    // Is useful in categories to construct cut strings out of each working point
    const std::map<lepID, std::string> map = { {L, "L"}, {T, "T"} };
  }
  
  // lepton Isolation
  namespace lepIso {
    enum lepIso{ no, L, T, Count };
    const std::array<lepIso, Count> it = {{ no, L, T }};
    const std::map<lepIso, std::string> map = { {no, "no"}, {L, "L"}, {T, "T"} };
  }

  // Combination of lepton ID + lepton Isolation for a single lepton
  uint16_t lepIDIso(const lepID::lepID& id, const lepIso::lepIso& iso);
  std::string lepIDIsoStr(const lepID::lepID& id, const lepIso::lepIso& iso);

  // Combination of lepton ID + lepton Isolation for a Dilepton object
  uint16_t leplepIDIso(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2);
  std::string leplepIDIsoStr(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2);

  // jet ID
  namespace jetID {
    enum jetID{ L, T, TLV, no, Count };
    const std::array<jetID, Count> it = {{ L, T, TLV, no }};
    const std::map<jetID, std::string> map = { {L, "L"}, {T, "T"}, {TLV, "TLV"}, {no, "no"} };
  }

  // B-tagging working points
  namespace btagWP {
    enum btagWP{ no, L, M, T, Count };
    const std::array<btagWP, Count> it = {{ no, L, M, T }};
    const std::map<btagWP, std::string> map = { {no, "no"}, {L, "L"}, {M, "M"}, {T, "T"} };
  }

  // jet combinatoric
  namespace jetPair {
    enum jetPair { ht, mh, pt, csv, jp, ptOverM, Count };
    const std::array<jetPair, Count> it = {{ ht, mh, pt, csv, jp, ptOverM }};
    const std::map<jetPair, std::string> map = { {ht, "ht"}, {mh, "mh"}, {pt, "pt"}, {csv, "csv"}, {jp, "jp"}, {ptOverM, "ptOverM"} }; 
  }

  // Combination of jet ID and B-tagging working point for one jet (not yet available) 
  uint16_t jetIDbtagWP(const jetID::jetID& id, const btagWP::btagWP& wp);
  std::string jetIDbtagWPStr(const jetID::jetID& id, const btagWP::btagWP& wp);
  
  // Combination of jet ID and B-tagging working points for two jets 
  uint16_t jetjetIDbtagWP(const jetID::jetID& id1, const jetID::jetID& id2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2);
  std::string jetjetIDbtagWPStr(const jetID::jetID& id1, const jetID::jetID& id2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2);

  // Combination of lepton ID, lepton Isolation, jet ID, B-tagging working points and jetPair ordering for a two-lepton-two-b-jets object
  uint16_t leplepIDIsojetjetIDbtagWPPair(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2, const jetID::jetID& jetid1, const jetID::jetID& jetid2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2, const jetPair::jetPair& jetpair);
  std::string leplepIDIsojetjetIDbtagWPPairStr(const lepID::lepID& id1, const lepIso::lepIso& iso1, const lepID::lepID& id2, const lepIso::lepIso& iso2, const jetID::jetID& jetid1, const jetID::jetID& jetid2, const btagWP::btagWP& wp1, const btagWP::btagWP& wp2, const jetPair::jetPair& jetpair);

}
