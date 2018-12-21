#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include <cstdlib>

#define RED 'r'
#define BLACK 'b'

using namespace std;

long pos1 , pos2 , pos3 , pos4 ;
string sval1 , sval2,sa;
long ival1,ival2,ia;
vector<string> line;
 vector<long>order;
 vector<long>ptree;
 vector<long>extree;
 vector<long>tottree;
 long tempid;
long temppos;
long temprem_time;
long temptot_time;
vector<string> jobid;
vector<string> arrayvalue;
int count=0;
int mark=0;
int mount=0;
long size;
void calculate();


struct node {
	long jobid;
	long extime;
	long total_time;
	struct node *left, *right, *p;
	char color;
};

typedef struct node *NODE;
struct node NIL;
NODE NULLPTR = &NIL;

NODE search(NODE root, long k)
{
	if (root == NULLPTR || root->jobid == k)
		return root;
	if (k >= root->jobid)
		return search(root->right, k);
	else
		return search(root->left, k);
}
void inorder(NODE n)//displays values in Red black tree
{
	if (n != NULLPTR) {
		inorder(n->left);
        cout<<n->jobid<<" ; ";
         cout<<n->extime<<" ; ";
          cout<<n->total_time<<" ; ";
		inorder(n->right);
	}
}
void porder(NODE n,long min,long max,long r1,long r2) //finds The jobid that exists in that range using BST
{


       if (n != NULLPTR) {
       if(r1<=min &&r2>=max)
       {

		porder(n->left,min,max,r1,r2);
		ptree.push_back(n->jobid);
		extree.push_back(n->extime);
		tottree.push_back(n->total_time);
		porder(n->right,min,max,r1,r2);
       }
       else if(r1>min&&r2>=max)
       {

           if(n->jobid<r1)
           {
            porder(n->right,min,max,r1,r2);
           }
           else if(n->jobid>r1 && n->jobid<r2)
           {
               porder(n->left,min,max,r1,r2);
               ptree.push_back(n->jobid);
               extree.push_back(n->extime);
		       tottree.push_back(n->total_time);
               porder(n->right,min,max,r1,r2);

           }
           else if(n->jobid==r1)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
             porder(n->right,min,max,r1,r2);
           }
           else if(n->jobid==r2)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
             porder(n->left,min,max,r1,r2);
           }
       }
       else if(r1<=min&&r2<=max)
       {

           if(n->jobid>r2)
           {
            porder(n->left,min,max,r1,r2);
           }
           else if(n->jobid>r1&& n->jobid<r2)
           {
               porder(n->left,min,max,r1,r2);
               ptree.push_back(n->jobid);
               extree.push_back(n->extime);
		       tottree.push_back(n->total_time);
               porder(n->right,min,max,r1,r2);

           }
           else if(n->jobid==r2)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
               porder(n->left,min,max,r1,r2);
           }
           else if(n->jobid==r1)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
               porder(n->right,min,max,r1,r2);
           }
       }
       else if(r1>min&&r2<max)
       {

          if(n->jobid>r1&& n->jobid<r2)
           {
               porder(n->left,min,max,r1,r2);
               ptree.push_back(n->jobid);
               extree.push_back(n->extime);
		       tottree.push_back(n->total_time);
               porder(n->right,min,max,r1,r2);

           }
           else if(n->jobid>r2)
           {
            porder(n->left,min,max,r1,r2);
           }
           else if(n->jobid<r1)
           {
            porder(n->right,min,max,r1,r2);
           }

           else if(n->jobid==r1)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
               porder(n->right,min,max,r1,r2);
           }
           else if(n->jobid==r2)
           {
             ptree.push_back(n->jobid);
             extree.push_back(n->extime);
		     tottree.push_back(n->total_time);
               porder(n->left,min,max,r1,r2);
           }
       }

	}
}
void rbupdateextime(NODE root, long jid, long time) //updates executed time of the job
{

    if (root->jobid == jid)
		root->extime=time;
	else if (jid >= root->jobid)
		 rbupdateextime(root->right,jid, time);
	else
		 rbupdateextime(root->left,jid, time);

}

void rbnextjob(NODE *root, long k)//finds if the next job exists
{
    ofstream createf;
    createf.open("output.txt",std::ios::app);
    long njob=-1;

    size_t i;

    NODE n =*root;
    while (n->right != NULLPTR)
		n = n->right;

    long max= n->jobid +300;

    for( i=0;i<order.size();i++){
            if(order.at(i)>k && order.at(i)<max)
            {
                njob=order.at(i);
                max=order.at(i);
            }

        }

       if(njob<0){
           createf<<"(0,0,0)";
           createf<<"\n";
        }
        else
        {

             NODE Z=search(*root,njob);
             createf<<"(";
             createf<<Z->jobid<<" ,";
             createf<<Z->extime<<", ";
             createf<<Z->total_time;
             createf<<")";
             createf<<"\n";

        }
createf.close();
}

void findjob(NODE n ,long k)// It searches the jobid like the BST and store it in a vector
{

       if (n != NULLPTR) {
            if(n->jobid==k)
            {
              order.push_back(n->jobid);
              if(n->left!=NULLPTR)
                order.push_back(n->left->jobid);
              if(n->right!=NULLPTR)
                order.push_back(n->right->jobid);


            }

       else if(n->jobid>k)
       {

           order.push_back(n->jobid);
           findjob(n->left,k);

       }
      else  if(n->jobid<k)
      {

           order.push_back(n->jobid);
           findjob(n->right,k);

      }

	}
}

void rbprevjob(NODE*root,long k )//checks if previous job exists
{

    ofstream createf;
    createf.open("output.txt",std::ios::app);
    long njob=-1;

    size_t i;

    long min=-1;
    for( i=0;i<order.size();i++){
       if(order.at(i)<k && order.at(i)>min)
            {
                njob=order.at(i);
                min=order.at(i);
            }

        }


       if(njob<0)
        {
           createf<<"(0,0,0)";
           createf<<"\n";

        }
        else{


             NODE Z=search(*root,njob);
             createf<<"(";
             createf<<Z->jobid<<" ,";
             createf<<Z->extime<<" , ";
             createf<<Z->total_time;
             createf<<")";
             createf<<"\n";

        }

createf.close();

}
NODE minimum(NODE root)//finds the minimum value in tree
{
	while (root->left != NULLPTR)
		root = root->left;
	return root;
}

NODE maxn(NODE root)//finds the maximum value in tree
{
	while (root->right != NULLPTR)
		root = root->right;
	return root;
}

void rbprintjob(NODE root, long k, int flag)//finds printjob(value)
{

    ofstream createf;
    createf.open("output.txt",std::ios::app);

   if(root!=NULLPTR)  {


	if ( root->jobid == k)
    {
        createf<<"(";
        createf<<root->jobid<<" , ";
        createf<<root->extime<<" , ";
        createf<<root->total_time;
        createf<<")";
        createf<<"\n";
        createf.close();
        flag=1;
        return;

    }

	else if (k >= root->jobid)
    {
     if(root->right!=NULLPTR)
       rbprintjob(root->right, k,flag);
      else
      {
         ofstream createf;
         createf.open("output.txt",std::ios::app);
         createf<<"(0,0,0)"<<endl;
         createf.close();
         return;
      }

    }

	else if (k < root->jobid)

		 {
           if(root->left!=NULLPTR)
             rbprintjob(root->left, k,flag);
             else
             {
                 ofstream createf;
                 createf.open("output.txt",std::ios::app);
                 createf<<"(0,0,0)"<<endl;
                 createf.close();
                 return;
             }
		 }
    return;
}

}



 long rbprintjob1(NODE root,long flag,long val1,long val2)// finds PrintJob(value,value)
{
    size_t i;
    size_t s;
    NODE x=NULL;
    NODE y=NULL;
    x = minimum(root);
    y = maxn(root);
    if(val2<x->jobid||val1>y->jobid)
    {

    }
    else
    {

        porder(root,x->jobid,y->jobid,val1,val2);
        ofstream createf;
        createf.open("output.txt",std::ios::app);
        s=ptree.size();

        for(i=0;i<s;i++)
        {
                flag=1;
                createf<<"(";
                createf<<ptree.at(i)<<" , ";
                createf<<extree.at(i)<<" , ";
                createf<<tottree.at(i);
                createf<<")";
                if(i<(s-1))
                createf<<",";
                if(i==(s-1))
                createf<<"\n";

        }

createf.close();
    }


return flag;

}


void lrotate(NODE *treeroot, NODE n)//left rotate the tree
 {
	NODE y = NULL;
	y= n->right;
	n->right = y->left;
	if (y->left != NULLPTR)
		y->left->p = n;
	y->p = n->p;
	if (n->p == NULLPTR)
		*treeroot = y;
	else if (n->p->left != n)
		n->p->right = y;
	else
		n->p->left = y;
	y->left = n;
	n->p = y;
}

void rrotate(NODE *treeroot, NODE n) {//right rotate the tree
	NODE x = NULL;
	x= n->left;
	n->left = x->right;
	if (x->right != NULLPTR)
		x->right->p = n;
	x->p = n->p;
	if (n->p == NULLPTR)
		*treeroot = x;
	else if (n->p->left != n)
		n->p->right = x;
	else
		n->p->left = x;
	x->right = n;
	n->p = x;
}

void rbinsertfix(NODE *treeroot, NODE z) {//rebalances the tree
	while (z->p->color == RED) {
		if (z->p == z->p->p->left) {
			NODE y = z->p->p->right;
            if (y->color != RED) {
               if (z == z->p->right) {
					z = z->p;
					lrotate(treeroot,z);
				}
				z->p->color = BLACK;
				z->p->p->color = RED;
				rrotate(treeroot,z->p->p);
            }
            else{
                z->p->color = BLACK;
				y->color = BLACK;
				z->p->p->color = RED;
				z = z->p->p;
			}


		}
		else {
			NODE y = z->p->p->left;
            if (y->color != RED) {
                    if (z == z->p->left) {
					   z = z->p;
					   rrotate(treeroot,z);
				    }
				    z->p->color = BLACK;
				    z->p->p->color = RED;
                    lrotate(treeroot,z->p->p);

            }
            else{
                z->p->color = BLACK;
				y->color = BLACK;
				z->p->p->color = RED;
				z = z->p->p;
            }


		}
	}
	(*treeroot)->color = BLACK;
}

void rbinsert(NODE *treeroot, long jid, long ttime) {//insert value in tree

	NODE Z = (NODE) malloc(sizeof(struct node));
	Z->jobid = jid;
	Z->extime=0;
	Z->total_time=ttime;
	NODE x = *treeroot;
	NODE y = NULLPTR;
	while (x != NULLPTR) {
		y = x;
		if (Z->jobid >= x->jobid)
			x = x->right;
		else
			x = x->left;
	}
	Z->p = y;
	if (y == NULLPTR)
		*treeroot = Z;
	else if (Z->jobid < y->jobid)
		y->left = Z;
	else
		y->right = Z;
	Z->left = NULLPTR;
	Z->right = NULLPTR;
	Z->color = RED;
	rbinsertfix(treeroot,Z);
}

void rbtransplant(NODE *treeroot, NODE u, NODE v) {
	if (u->p == NULLPTR)
		*treeroot = v;
	else if (u == u->p->left)
		u->p->left = v;
	else
		u->p->right = v;
	v->p = u->p;
}

void rbdeletefix(NODE *treeroot, NODE x) {//rebalances after delete
	while (x != *treeroot && x->color == BLACK) {
		if (x == x->p->left) {
			NODE w = x->p->right;
			if (w->color == RED) {
				w->color = BLACK;
				x->p->color = RED;
				lrotate(treeroot,x->p);
				w = x->p->right;
			}
			if (w->left->color == BLACK && w->right->color == BLACK) {
				w->color = RED;
				x = x->p;
			}
			else {
			 	if (w->right->color == BLACK) {
					w->left->color = BLACK;
					w->color = RED;
					rrotate(treeroot,w);
					w = x->p->right;
				}
				w->color = x->p->color;
				x->p->color = BLACK;
				w->right->color = BLACK;
				lrotate(treeroot,x->p);
				x = *treeroot;
			}
		}
		else {
			NODE w = x->p->left;
			if (w->color == RED) {
				w->color = BLACK;
				x->p->color = RED;
				rrotate(treeroot,x->p);
				w = x->p->left;
			}
			if (w->left->color == BLACK && w->right->color == BLACK) {
				w->color = RED;
				x = x->p;
			}
			else {
				if (w->left->color == BLACK) {
					w->right->color = BLACK;
					w->color = RED;
					lrotate(treeroot,w);
					w = x->p->left;
				}
				w->color = x->p->color;
				x->p->color = BLACK;
				w->left->color = BLACK;
				rrotate(treeroot,x->p);
				x = *treeroot;
			}
		}
	}
	x->color = BLACK;
}

void rbdelete(NODE *treeroot, long z) {//deletes value in tree
	NODE Z = search(*treeroot, z);
	if (Z == NULLPTR) {

		return;
	}
	NODE y = NULL;
	y=Z;
	long yoc = y->color;
	NODE x;
	if (Z->left == NULLPTR) {
		x = Z->right;
		rbtransplant(treeroot,Z,Z->right);
	}
	else if (Z->right == NULLPTR) {
		x = Z->left;
		rbtransplant(treeroot,Z,Z->left);
	}
	else {
		y = minimum(Z->right);
		yoc = y->color;
		x = y->right;
		if (y->p == Z)
			x->p = y;
		else {
			rbtransplant(treeroot,y,y->right);
			y->right = Z->right;
			y->right->p = y;
		}
		rbtransplant(treeroot,Z,y);
		y->left = Z->left;
		y->left->p = y;
		y->color = Z->color;
	}
	if (yoc == BLACK)
		rbdeletefix(treeroot,x);
}


template <class T>
class MinHeap {
  long getLeftChild(long parent);
  long getRightChild(long parent);

public:
  MinHeap();
  vector<T> pos;
  vector<T> jid;
  vector<T>rem_time;
  vector<T>total_time;
  vector<long>fin_jid;
  vector<long>fin_time;
  void heapify();
  void inserth(T ,T );
  void swap(long child, long parent);
  T remove();
  long getSize();
  void setSize(long x);
  void afterdel();
  void display();

};

template <class T>
MinHeap<T> :: MinHeap(){

}

template <class T>
long MinHeap<T> :: getSize(){//get size of heap
  return pos.size();
}

template <class T>
void MinHeap<T> :: setSize(long x){//resize the heap

  pos.resize(x);

}

template <class T>
void MinHeap<T> :: display(){//display min heap
size_t i=0;
size_t v=0;
 cout<<"display--"<<"\n";
for(i=0;i<jid.size();i++)
 {
   cout<<jid[i]<<"  "<<pos[i]<<"\n";
 }
 cout<<"complete "<<endl;
 if(fin_jid.size()<=(size_t)0)
 {
     cout<<"NULL"<<endl;
 }
 else{
    for( v=0;v<fin_jid.size();v++)
     {
      cout<<fin_jid[v]<<"  "<<fin_time[v]<<"\n";
    }
 }
}



template <class T>
void MinHeap<T>::swap(long child, long parent) {//swap values in min heap
  T temp;
  T temp1;
  T temp3;
  T  temp4;
  temp = pos[child];
  temp1=jid[child];
  temp3= total_time[child];
  temp4= rem_time[child];
  pos[child] = pos[parent];
  jid[child] = jid[parent];
  total_time[child]=total_time[parent];
  rem_time[child]=rem_time[parent];
  pos[parent] = temp;
  jid[parent] = temp1;
  total_time[parent]=temp3;
  rem_time[parent]=temp4;
}
template <class T>
long MinHeap<T> :: getLeftChild(long parent){//gets left child
  return 2*parent +1;
}
template <class T>
long MinHeap<T> :: getRightChild(long parent){//gets right child
  return 2 * parent + 2;
}
template <class T>
void MinHeap<T> :: inserth(T val,T val1) {// insert values in min heap
int mount=1;
size_t i=jid.size();
pos[i]=0;
total_time.push_back(val);
rem_time.push_back(val);
jid.push_back(val1);
if(mark==0)
  heapify();
else if(mark ==1)
 {
    size_t l;
    long temp1,temp2,temp3,temp4;
    for(  l=0;l< (jid.size()-1);l++)
         {
               jid[l]=jid[l+1];
               //pos[l]=pos[l+1];
               rem_time[l]=rem_time[l+1];
               total_time[l]=total_time[l+1];

           }
    for(  l=0;l< (pos.size()-1);l++)
         pos[l]=pos[l+1];

    jid.pop_back();
    rem_time.pop_back();
    total_time.pop_back();
    afterdel();
    heapify();
    jid.push_back(tempid);
    rem_time.push_back(temprem_time);
    total_time.push_back(temptot_time);

    pos[i]=temppos;
    swap(0,i);

   }

}
template <class T>
void MinHeap <T>:: heapify() {//heapify the min heap

  long child=0;
  size_t i=0;

   i=jid.size();
   child=i-1;
  long parent= 0;

  if (child % 2 == 0)
	parent= (child /2 ) -1;
  else
	parent= child/2;


  while (pos[child] < pos[parent] && child >=0 && parent >= 0) {

	swap(child, parent);
	child = parent;
	if (child % 2 == 0)
	parent= (child /2 ) -1;
   else
	parent= child/2;
  }

}


template <class T>
T MinHeap<T> :: remove() {//remove element

  long child=0;
  size_t i;
  size_t j;


   i=jid.size();
  child=i-1;
  swap(child, 0);

  T value = pos[child];


  for( j=i-1;j<pos.size()-1;j++)
    pos[j]=pos[j+1];

  for( j=i-1;j<jid.size()-1;j++)
    jid[j]=jid[j+1];
  for( j=i-1;j<jid.size()-1;j++)
    total_time[j]=total_time[j+1];
  for( j=i-1;j<jid.size()-1;j++)
    rem_time[j]=rem_time[j+1];

  pos.pop_back();
  jid.pop_back();
  total_time.pop_back();
  rem_time.pop_back();
  afterdel();

  return value;


}


template <class T>
void MinHeap<T> :: afterdel() {//rebalance after delete
  long parent = 0;

  while (1) {

	long left = getLeftChild(parent);
	long right = getRightChild(parent);


    long length = jid.size();
	long largest = parent;

	if (left < length && pos[left] < pos[largest])
	  largest = left;

	if (right < length && pos[right] < pos[largest])
	  largest = right;

	if (largest != parent) {
	  swap(largest, parent);
	  parent = largest;
	}
	else
	  break;


  }

}



void calculate()
{

    long size1=0;
    MinHeap<long> heap;

    NIL.left = NIL.right = NIL.p = NULLPTR;
	NIL.color = BLACK;
	NODE tree = NULLPTR;


    static long insertcount=0;
    static long progcounter=0;
    size=line.size();

    for (long i=0;i<size;i++)
    {

        string func;
        long pos1=pos2=0;
        pos1=line[i].find(":");
        pos2=line[i].find("(");
        func =line[i].substr(pos1+2,pos2-(pos1+2));
        sa=line[i].substr(0,pos1);
        ia=atoi(sa.c_str());

        while(progcounter<ia)
                {
                    if(count==0 && mark==0 && (tree!=NULLPTR))
                      {
                         tempid=heap.jid[0];
                         temppos=heap.pos[0];
                         temprem_time=heap.rem_time[0];
                         temptot_time=heap.total_time[0];
                         count= temprem_time;
                         if(temprem_time>5)
                            count=5;
                         mark=1;
                       }
                    if(count>0)
                      {
                         heap.pos[0]=heap.pos[0]+1;
                         temppos=temppos+1;
                         rbupdateextime(tree,heap.jid[0],heap.pos[0]);
                         heap.rem_time[0]=heap.rem_time[0]-1;
                         temprem_time=temprem_time-1;
                         count=count-1;
                      }
                   if(count==0 && mark==1)
                     {
                        if(heap.rem_time[0]==0)
                           {

                             rbdelete(&tree, heap.jid[0]);
                             heap.remove();
                             mark =0;
                             count=0;

                            }
                   else{
                       rbupdateextime(tree,heap.jid[0],heap.pos[0]);
                       size_t i;
                       i=heap.jid.size();
                       if(mount ==1)
                        {
                      heap.swap(0,i-1);
                      mount=0;

                       }

                  else{
                        size_t l;
                         for(  l=0;l< (heap.jid.size()-1);l++)
                           {
                             heap.jid[l]=heap.jid[l+1];

                             heap.rem_time[l]=heap.rem_time[l+1];
                              heap.total_time[l]= heap.total_time[l+1];

                             }
                         for( l=0;l< (heap.pos.size()-1);l++)
                              heap.pos[l]= heap.pos[l+1];

                         heap.jid.pop_back();
                         heap.rem_time.pop_back();
                         heap.total_time.pop_back();
                         heap.afterdel();
                         heap.heapify();
                         heap.jid.push_back(tempid);
                         heap.rem_time.push_back(temprem_time);
                         heap.total_time.push_back(temptot_time);
                         i=heap.jid.size();
                         heap.pos[i-1]=temppos;

                  }

                heap.afterdel();
                heap.heapify();
                mark=0;
                count=0;
                 }

            }

           progcounter=progcounter+1;
      }

        if (progcounter==ia)
        {

            if((func.compare("Insert"))==0){

            insertcount++;
            if(size1<insertcount)
            {
                if(size1==0)
                    size1=1;
                else
                    size1=size1*2;
                 heap.setSize(size1);

            }
            pos3=pos4=0;
            ival1=ival2=0;
            pos3=line[i].find(",");
            pos4=line[i].find(")");

            sval1= line[i].substr(pos3+1,pos4-(pos3+1));
            ival1= atoi(sval1.c_str());

            sval2= line[i].substr(pos2+1,pos3-(pos2+1));
            ival2= atoi(sval2.c_str());

            arrayvalue.push_back(sval1);
            jobid.push_back(sval2);

             heap.inserth(ival1,ival2); //heap insert job id value

             rbinsert(&tree,ival2,ival1);// red black tree insert  job id arrival time


            progcounter=progcounter+1;

           }
           if((func.compare("PrintJob"))==0)
           {

              long pos3=pos4=0;
              pos3=line[i].find(",");
              pos4=line[i].find(")");

              if(pos3<=0){
                 string valp=line[i].substr(pos2+1,pos4-(pos2+1));//checks if it is a print job with Single value or double
                 long ivalp=atoi(valp.c_str());

                 rbprintjob(tree,ivalp,0);//PrintJob(val)

                }
              else{//printJob(val1, val2)

                  string valp1= line[i].substr(pos2+1,pos3-(pos2+1));
                  long ivalp1 = atoi(valp1.c_str());
                  string valp2= line[i].substr(pos3+1,pos4-(pos3+1));
                  long ivalp2 = atoi(valp2.c_str());

                  long flag=0;
                  ptree.clear();
                  extree.clear();
                  tottree.clear();

                   flag= rbprintjob1(tree,flag,ivalp1,ivalp2);

                if (flag==0)
                 {
                     ofstream createf;
                     createf.open("output.txt",std::ios::app);
                     createf<<"(0,0,0)"<<endl;
                     createf.close();

                 }

              }

                progcounter=progcounter+1;
           }
            if((func.compare("NextJob"))==0)
           {

               long pos4=0;
               pos4=line[i].find(")");
               string valp=line[i].substr(pos2+1,pos4-(pos2+1));
               long ivalp=atoi(valp.c_str());
               order.clear();
               findjob(tree,ivalp);
               rbnextjob(&tree,ivalp);
               progcounter=progcounter+1;
           }

           if((func.compare("PreviousJob"))==0)
           {

               long pos4=0;
               pos4=line[i].find(")");
               string valp=line[i].substr(pos2+1,pos4-(pos2+1));
               long ivalp=atoi(valp.c_str());
               order.clear();
               findjob(tree,ivalp);
               rbprevjob(&tree,ivalp);
               progcounter=progcounter+1;
           }
                  //Dispatching job an after an operation
                    if(count==0 && mark==0 && (tree!=NULLPTR))
                      {
                       tempid=heap.jid[0];
                       temppos=heap.pos[0];
                       temprem_time=heap.rem_time[0];
                       temptot_time=heap.total_time[0];
                       count= temprem_time;
                       if(temprem_time>5)
                        count=5;
                       mark=1;
                       }
                   if(count>0)
                      {
                         heap.pos[0]=heap.pos[0]+1;
                         temppos=temppos+1;
                         rbupdateextime(tree,heap.jid[0],heap.pos[0]);
                         heap.rem_time[0]=heap.rem_time[0]-1;
                         temprem_time=temprem_time-1;
                         count=count-1;
                      }
                  if(count==0 && mark==1)
                     {
                        if(heap.rem_time[0]==0)
                           {

                           rbdelete(&tree, heap.jid[0]);
                           heap.remove();
                           mark =0;
                           count=0;
                 }
                 else{
                  rbupdateextime(tree,heap.jid[0],heap.pos[0]);
                  size_t i;
                  i=heap.jid.size();
                  if(mount ==1)
                  {
                      heap.swap(0,i-1);
                      mount=0;

                  }


                  else{
                        size_t l;
                         for(  l=0;l< (heap.jid.size()-1);l++)
                           {
                             heap.jid[l]=heap.jid[l+1];
                             heap.rem_time[l]=heap.rem_time[l+1];
                             heap.total_time[l]= heap.total_time[l+1];

                             }
                       for(  l=0;l< (heap.pos.size()-1);l++)
                            heap.pos[l]= heap.pos[l+1];

                       heap.jid.pop_back();
                       heap.rem_time.pop_back();
                       heap.total_time.pop_back();
                       heap.afterdel();
                       heap.heapify();
                       heap.jid.push_back(tempid);
                       heap.rem_time.push_back(temprem_time);
                       heap.total_time.push_back(temptot_time);
                       i=heap.jid.size();
                       heap.pos[i-1]=temppos;


                  }

                 heap.afterdel();
                 heap.heapify();
                 mark=0;
                 count=0;
            }

        }


    }


}


// file is completely read but there is still values in heap

    long y=1;
    while(y==1){

              if(count==0 && mark==0 && (tree!=NULLPTR))
                      {
                       tempid=heap.jid[0];
                       temppos=heap.pos[0];
                       temprem_time=heap.rem_time[0];
                       temptot_time=heap.total_time[0];
                       count= heap.rem_time[0];
                       if(heap.rem_time[0]>5)
                        count=5;
                       mark=1;
                       }
                  if(count>0)
                   {
                    heap.pos[0]=heap.pos[0]+1;
                    temppos=temppos+1;
                    rbupdateextime(tree,heap.jid[0],heap.pos[0]);
                    heap.rem_time[0]=heap.rem_time[0]-1;
                    temprem_time=temprem_time-1;
                    count=count-1;
                    }
               if(count==0 && mark==1)
                {
                  if(heap.rem_time[0]==0)
                   {
                 rbdelete(&tree, heap.jid[0]);
                 heap.remove();

                  mark =0;
                  count=0;

                 }
                 else{
                  rbupdateextime(tree,heap.jid[0],heap.pos[0]);

                  size_t i;
                  i=heap.jid.size();
                  if(mount ==1)
                  {
                      heap.swap(0,i-1);
                      mount=0;

                  }


                  else{
                        size_t l;
                         for(  l=0;l< (heap.jid.size()-1);l++)
                           {
                             heap.jid[l]=heap.jid[l+1];
                             heap.rem_time[l]=heap.rem_time[l+1];
                             heap.total_time[l]= heap.total_time[l+1];

                             }
                       for(  l=0;l< (heap.pos.size()-1);l++)
                            heap.pos[l]= heap.pos[l+1];

                       heap.jid.pop_back();
                       heap.rem_time.pop_back();
                       heap.total_time.pop_back();
                       heap.afterdel();
                       heap.heapify();
                       heap.jid.push_back(tempid);
                       heap.rem_time.push_back(temprem_time);
                       heap.total_time.push_back(temptot_time);
                       i=heap.jid.size();
                       heap.pos[i-1]=temppos;
                  }

                 heap.afterdel();
                 heap.heapify();
                 mark=0;
                count=0;
                 }

                }

    if(heap.jid.size()==0)//If there is no value in heap then the while loop terminates
       y=0;

    }
    }


int main(int argc,char *argv[]){

    string data;
    string val;
    string fname;

     ofstream createf;
    createf.open("output_file.txt", std::ios::out | std::ios::trunc);//If a file name output exist remove it and create a new one.
    createf.close();

    fname=argv[1];
    ifstream readfile;
    fname="./"+fname;
	readfile.open(fname.c_str());
	while(!readfile.eof())
    {

      getline(readfile,data);
      line.push_back(data);
    }

    calculate();


   return 0;

}

