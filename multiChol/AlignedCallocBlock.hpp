#pragma once

#include <cstdlib>
#include <new> // For std::bad_alloc

template <typename T>
class AlignedCallocBlock {
public:
   AlignedCallocBlock(AlignedCallocBlock const&) = delete;
   void operator=(AlignedCallocBlock const&) = delete;
   AlignedCallocBlock(size_t nmemb) {
      nmemb += 32/sizeof(T); // Ensure space to align to 32-byte boundary
      base_ = calloc(nmemb, sizeof(T)); // Allocate memory
      if(!base_) throw std::bad_alloc();
      unsigned long cbase = (unsigned long) base_; // Cast to char for pointer arithmetic
      cbase = 32*(cbase / 32) + 32; // Must be 32-byte aligned and within blk
      ptr_ = (T *) cbase; // Cast to desired type
   }
   ~AlignedCallocBlock() {
      free(base_);
   }
   T *get_ptr() const {
      return ptr_;
   }
private:
   void *base_;
   T *ptr_;
};

