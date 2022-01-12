#include <InitDictRandomly.h>
#include <Shuffle.h>

#include <random>

using namespace ttk;

void InitRandomDict::execute(std::vector<Diagram> &DictDiagrams,
                             const std::vector<Diagram> &datas,
                             const int nbAtom,
                             const int seed) {
  int nDiags = datas.size();
  DictDiagrams.resize(nbAtom);
  std::vector<int> indices(nDiags);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 random_engine{};
  random_engine.seed(seed);
  ttk::shuffle(indices, random_engine);
  for(int i = 0; i < nbAtom; ++i) {
    const Diagram &atom = datas[indices[i]];
    DictDiagrams[i] = atom;
  }
}
