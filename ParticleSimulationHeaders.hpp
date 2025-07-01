/*----------
  ver 250630
  ----------*/

#ifndef ParticleSimulationHeaders_hpp
#define ParticleSimulationHeaders_hpp

#ifndef vs

#include<bits/stdc++.h>

#else

#define _SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING

#include <cassert>

#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdalign>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>

#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#endif

//#include"numcxx_complex.hpp"

#include"./BasicTools.hpp"
#include"./AgentRealization.hpp"
#include"./StateProcessing.hpp"
#include"./CollectiveMotionRealization.hpp"
#include"./InteractionRealization.hpp"
#include"./DataProcessing.hpp"
#include"./ClusterUnionFind.hpp"

#include"./ReducedModel.hpp"

#include"./ActiveIsingModel.hpp"

#include"./VicsekModel.hpp"
#include"./ContinuousVicsekModel.hpp"
#include"./NematicAlgnPP.hpp"
#include"./ActiveNematics.hpp"

#ifdef CheckingHeader
HeaderChecker(ParticleSimulationHeaders_hpp);
#endif //CheckingHeader

#endif //Particlesimulationheaders
