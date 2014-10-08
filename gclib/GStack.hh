//basic template for a Stack of object pointers (implemented as a linked list)
#include "GBase.h"

template <class OBJ> class GStack {
 protected:
   struct StackOBJ {
      OBJ* obj;
      StackOBJ* prev;
      };
   int fCount; //total number of elements in stack
   StackOBJ* base;
   StackOBJ* top;
 public:
   GStack(OBJ* po=NULL) {
      base=NULL;
      top=NULL;
      fCount=0;
      if (po!=NULL) Push(po);
      }
   ~GStack() {
      while (fCount>0) Pop();
      }
   bool isEmpty() { return fCount==0; }
   int Size() { return fCount; }
   int Count() { return fCount; }
   OBJ* Pop() {
      if (top==NULL) return NULL;
      fCount--;
      StackOBJ* ctop=top;
      if (top==base) base=NULL;
      OBJ* r=top->obj;
      top=top->prev;
      GFREE(ctop);
      return r;
      }
   OBJ* Push(OBJ* o) {
      fCount++;
      StackOBJ* ctop=top; //could be NULL
      GMALLOC(top, sizeof(StackOBJ));
      top->obj=o;
      top->prev=ctop;
      if (base==NULL) base=top;
      return o;
      }
   OBJ* Top() { return ((top==NULL)? NULL : top->obj); }
   OBJ* Base() { return ((base==NULL)? NULL : base->obj); }
};
