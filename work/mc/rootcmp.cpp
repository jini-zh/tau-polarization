#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <cmath>
#include <cstddef>
#include <cstring>

#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TLeafD.h>
#include <TTree.h>

class UnknownType: public std::exception {
  public:
    UnknownType(const char* type);
    const char* what() const noexcept;
  private:
    std::string message_;
};

class bail {
  public:
    operator bool();
    template <typename T> bail& operator<<(const T&);
};

UnknownType::UnknownType(const char* type):
  message_(std::string("Unknown ROOT type: ") + type) {};

const char*
UnknownType::what() const noexcept {
  return message_.c_str();
};

bail::operator bool() {
  return false;
};

template <typename T> bail& bail::operator<<(const T& x) {
  std::cout << x;
  return *this;
};

size_t leaves_memory_size(TBranch& branch) {
  TObjArray* leaves = branch.GetListOfLeaves();
  if (!leaves) return 0;

  size_t result = 0;
  for (TObject* leaf: *leaves)
    result += dynamic_cast<TLeaf*>(leaf)->GetLenType();

  return result;
};

size_t branches_memory_size(TObjArray* branches) {
  if (!branches) return 0;

  size_t result = 0;
  for (TObject* b: *branches) {
    TBranch& branch = dynamic_cast<TBranch&>(*b);
    result += leaves_memory_size(branch);
    result += branches_memory_size(branch.GetListOfBranches());
  };

  return result;
};

unsigned char*
set_branches_addresses(TObjArray* branches, unsigned char* buffer) {
  if (!branches) return buffer;

  for (TObject* b: *branches) {
    TBranch& branch = dynamic_cast<TBranch&>(*b);
    branch.SetAddress(buffer);
    buffer += leaves_memory_size(branch);
    buffer = set_branches_addresses(branch.GetListOfBranches(), buffer);
  };

  return buffer;
};


// Compares the contents of two objects
bool equal(TObject& a, TObject& b);

// Compares the structure of two objects.
bool equal_structure(TObject& a, TObject& b);

// Compares the contents of two objects assuming the same structure.
bool equal_data(TObject& a, TObject& b);

bool equal(TKey& a, TKey& b) {
  if (std::strcmp(a.GetName(), b.GetName()) != 0)
    return bail() << "Different keys: `"
      << a.GetName() << "' and `" << b.GetName() << "'\n";
  if (!equal(*a.ReadObj(), *b.ReadObj()))
    return bail() << "  in " << a.GetName() << '\n';
};

bool equal(TList& a, TList& b) {
  TObjLink* la = a.FirstLink();
  TObjLink* lb = b.FirstLink();
  for (int i = 0;; ++i) {
    if (la)
      if (lb) { // la && lb
        if (!equal(*la->GetObject(), *lb->GetObject()))
          return bail() << "  at position " << i << '\n';
      } else // la && !lb
        return bail() << "First list has less elements (" << i << ") than the second list\n";
      else if (lb) // !la && lb
        return bail() << "First list has more elements (" << i << ") than the second list\n";
      else // !la && !lb
        return true;
    la = la->Next();
    lb = lb->Next();
  };
  return true;
};

bool equal_structure(TObjArray& a, TObjArray& b) {
  Long64_t n  = a.GetEntries();
  Long64_t nb = b.GetEntries();
  if (n != nb)
    return bail() << "Different array sizes: " << n << " != " << nb << '\n';

  for (Long64_t i = 0; i < n; ++i)
    if (!equal_structure(*a[i], *b[i]))
      return bail() << "  at position " << i << '\n';

  return true;
};

bool equal_data(TObjArray& a, TObjArray& b) {
  for (Long64_t i = 0; i < a.GetEntries(); ++i)
    if (!equal_data(*a[i], *b[i]))
      return bail() << "  at position " << i << '\n';
  return true;
};

bool equal(TH1& a, TH1& b) {
  if (std::strcmp(a.GetName(), b.GetName()) != 0)
    return bail() << "Different histograms: `"
      << a.GetName() << "' and `" << b.GetName() << "'\n";

  if (a.GetBufferSize() != b.GetBufferSize())
    return bail() << "Histograms " << a.GetName() << " are of different size: "
      << a.GetBufferSize() << " != " << b.GetBufferSize() << '\n';

  if (std::memcmp(a.GetBuffer(), b.GetBuffer(), a.GetBufferSize()) != 0)
    return bail() << "Different contents of histogram " << a.GetName() << '\n';

  return true;
};

bool equal_structure(TLeaf& a, TLeaf& b) {
  if (std::strcmp(a.GetName(), b.GetName()) != 0)
    return bail() << "Different leaves: `"
      << a.GetName() << "' and `" << b.GetName() << "'\n";

  if (std::strcmp(a.GetTypeName(), b.GetTypeName()) != 0)
    return bail() << "Different leaf types: `"
      << a.GetTypeName() << "' and `" << b.GetTypeName()
      << "\n  in leaf " << a.GetName() << '\n';

  return true;
};

bool equal_data(TLeafD& a, TLeafD& b) {
  if (a.GetLen() != b.GetLen())
    return bail() << "Different arrays sizes in leaf " << a.GetName() << ": "
      << a.GetLen() << " != " << b.GetLen() << '\n';

  if (std::memcmp(a.GetValuePointer(), b.GetValuePointer(), a.GetLen() * a.GetLenType()) != 0)
    for (Int_t i = 0; i < a.GetLen(); ++i)
      if (a.GetValue(i) != b.GetValue(i)) {
        double x = a.GetValue(i);
        double y = b.GetValue(i);
        double d = std::abs(x / y - 1);
//        if (x != y && x != 0 && d > 10 * std::numeric_limits<double>::epsilon())
        if (x != y && x != 0 && d > 1e-10)
        {
          bail log;
          log << x << " != " << y;
          if (d < 1e-7 || true) log << " (" << d << ")";
          if (a.GetLen() > 0) log << " at position " << i;
          log << "\n  in leaf " << a.GetName() << '\n';
          return false;
        };
      };

  return true;
};

bool equal_data(TLeaf& a, TLeaf& b) {
  if (a.GetLen() != b.GetLen())
    return bail() << "Different arrays sizes in leaf " << a.GetName() << ": "
      << a.GetLen() << " != " << b.GetLen() << '\n';

  if (std::memcmp(a.GetValuePointer(), b.GetValuePointer(), a.GetLen() * a.GetLenType()) != 0)
    for (Int_t i = 0; i < a.GetLen(); ++i)
      if (a.GetValue(i) != b.GetValue(i)) {
        bail log;
        log << a.GetValue(i) << " != " << b.GetValue(i);
        if (a.GetLen() > 0) log << " at position " << i;
        log << "\n  in leaf " << a.GetName() << " (" << a.ClassName() << ")\n";
        return false;
      };

  return true;
};

bool equal_structure(TBranch& a, TBranch& b) {
  if (std::strcmp(a.GetName(), b.GetName()) != 0)
    return bail() << "Different branches: `"
      << a.GetName() << "' and `" << b.GetName() << "'\n";

  if (!equal_structure(*a.GetListOfLeaves(), *b.GetListOfLeaves()))
    return bail() << "  in branch " << a.GetName() << '\n';

  if (!equal_structure(*a.GetListOfBranches(), *b.GetListOfBranches()))
    return bail() << "  in branch " << a.GetName() << '\n';

  return true;
};

bool equal_data(TBranch& a, TBranch& b) {
  if (!equal_data(*a.GetListOfLeaves(), *b.GetListOfLeaves()))
    return bail() << "  in branch " << a.GetName() << '\n';

  if (!equal_data(*a.GetListOfBranches(), *b.GetListOfBranches()))
    return bail() << "  in branch " << a.GetName() << '\n';

  return true;
};

bool equal(TTree& a, TTree& b) {
  if (std::strcmp(a.GetName(), b.GetName()) != 0)
    return bail() << "Different trees: `"
      << a.GetName() << "' and `" << b.GetName() << "'\n";

  if (!equal_structure(*a.GetListOfBranches(), *b.GetListOfBranches()))
    return bail() << "  when comparing structures of tree " << a.GetName() << '\n';

  if (a.GetEntries() != b.GetEntries())
    return bail() << "Different number of entries in tree " << a.GetName() << ": "
      << a.GetEntries() << " != " << b.GetEntries() << '\n';
  if (a.GetEntries() == 0) return true;

  for (Long64_t i = 0; i < a.GetEntries(); ++i) {
    a.GetEntry(i);
    b.GetEntry(i);
    if (!equal_data(*a.GetListOfBranches(), *b.GetListOfBranches()))
      return bail() << "  in entry " << i << " in tree " << a.GetName() << '\n';
  };

  return true;
};

bool equal(TDirectoryFile& a, TDirectoryFile& b) {
  return equal(*a.GetListOfKeys(), *b.GetListOfKeys());
};

// The classes below should go in reverse inheritance order:
// from children to parents
#define apply_macro_to_root_classes(macro) \
  macro(TDirectoryFile) \
  macro(TKey) \
  macro(TTree) \
  macro(TBranch) \
  macro(TLeafD) \
  macro(TLeaf) \
  macro(TH1)

#define dispatch_f(type, function) \
  { \
    type* pa = dynamic_cast<type*>(&a); \
    if (pa) return function(*pa, dynamic_cast<type&>(b)); \
  };

bool equal(TObject& a, TObject& b) {
  if (std::strcmp(a.ClassName(), b.ClassName()) != 0)
    return bail() << "Different object types: `" << a.ClassName() << "' and `"
      << b.ClassName() << "'\n";

#define dispatch(type) dispatch_f(type, equal)
  apply_macro_to_root_classes(dispatch);
#undef dispatch

  throw UnknownType(a.ClassName());
};

bool equal_structure(TObject& a, TObject& b) {
  if (std::strcmp(a.ClassName(), b.ClassName()) != 0)
    return bail() << "Different data types: `"
      << a.ClassName() << "' and `" << b.ClassName() << "'\n";

#define dispatch(type) dispatch_f(type, equal_structure)
  apply_macro_to_root_classes(dispatch);
#undef dispatch

  throw UnknownType(a.ClassName());
};

bool equal_data(TObject& a, TObject& b) {
  if (std::strcmp(a.ClassName(), b.ClassName()) != 0)
    return bail() << "Different data types: `"
      << a.ClassName() << "' and `" << b.ClassName() << "'\n";

#define dispatch(type) dispatch_f(type, equal_data)
  apply_macro_to_root_classes(dispatch);
#undef dispatch

  throw UnknownType(a.ClassName());
};

#undef dispatch_f
#undef apply_macro_to_root_classes

//template <typename T>
//typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
//equal(T a, T b) {
//  if (a != b)
//    return bail() << a << " != " << b << '\n';
//  return true;
//};

#if 0
template <typename T>
typename std::enable_if<!std::is_pointer<T>::value, bool>::type
equal(T& a, T& b) {
  if (std::memcmp(&a, &b, sizeof(T)) != 0)
    return bail() << "There is a difference in data\n";
  return true;
};
#endif

int main(int argc, char** argv) {
  TFile a(argv[1]);
  TFile b(argv[2]);
  return !equal(a, b);
};
