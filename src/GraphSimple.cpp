//============================================================================
// Name        : GraphSimple.cpp
// Author      : stam
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <ctime>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <limits>
#include <queue>
#include <array>
#include <set>
#include<cstdio>

using namespace std;


class Edge{

public:
	int firstNode;
	int secondNode;
	double weight;
	bool visited;

public:
	Edge(int from,int to,double cost)
	{
		firstNode=from;
		secondNode=to;
		weight=cost;
		visited=false;
	}

	int getFirstNode()
	{
		return firstNode;
	}

	int getSecondNode()
	{
		return secondNode;
	}

	double getWeight()
	{
		return weight;
	}
	void setWeight(double thevalue)
	{
		weight=thevalue;
	}
	void setVisited(bool value)
	{
		visited=value;
	}
	bool getVisited()
	{
		return visited;
	}
	 bool operator<(const  Edge& second) const {
	    return firstNode+secondNode < second.firstNode +second.secondNode ;
	  }
	 bool operator==(const  Edge& second) const {
	 	    return (firstNode==second.secondNode )&&(secondNode == second.secondNode) ;
	 	  }

};

class Node
{
private:
	int id;
	std::vector<Edge> adjNodeslist;
public:
	Node(int nodeId)
	{
		id=nodeId;
	}
	Node()
	{

	}
	~Node()
	{
		adjNodeslist.clear();
	}
	inline int getId()
	const
	{
		return id;
	}
	int addAdjEdge(Edge& a){

		int found=0;
		for (int i=0;i<adjNodeslist.size();i++)
		{
			if  ( ( a.getFirstNode()==adjNodeslist.at(i).getFirstNode() )&& (a.getSecondNode()==adjNodeslist.at(i).getSecondNode() ) )
				found=1;
		}
		if(!found){
			adjNodeslist.push_back(a);
			return 1;
		}else{
			return 0;
		}

	}

	std::vector<Edge>& getAdjEdgeList()
	{
		return adjNodeslist;
	}
	void printAdjList()
	{
		for (int i=0 ; i < adjNodeslist.size() ; i++ )
		{
			Edge edge = adjNodeslist[i];
			cout << id << "->" << edge.getSecondNode();
		}
	}
	bool operator<( const Node& rhs ) const
	        { return  id < rhs.getId() ;}
};

class distanceToNode{
private:
	Node current;
	int predecessor;
	bool known;
public:
	double distance;
	distanceToNode()
	{
		known=false;
		distance=numeric_limits<double>::max();
	}
	double getDistance()
	const
	{
		return distance;
	}
	void inline setDistance(double new_distance)
	{
		distance=new_distance;
	}
	void addToDistance(double increment)
	{
		distance+=increment;
	}
	void setPredecessor(int aNode)
	{
		predecessor=aNode;
	}
	int getPredecessor()
	{
		return predecessor;
	}
	void addCurrent(Node& aNode)
	{
		current=aNode;
	}
	Node getCurrent()
	{
		return current;
	}
};

struct CompdistanceToNode
{
  bool operator()(const  distanceToNode first, const  distanceToNode second) const {
    return first.distance > second.distance ;
  }
};

class Graph {
private:
	int numberOfNodes;
	vector<Node> listOfNodes;
public:
	//constructor
	Graph()
	{
		numberOfNodes=0;
	}
	Graph(int noOfNodes)
	{
		numberOfNodes=noOfNodes;
		for(int i=0 ;i < noOfNodes ;i++)
		{
			Node aNode(i);
			listOfNodes.push_back(aNode);
		}

	}
	//destructor
	// ~Graph();
	inline vector<Node>&  getListOfNodes()
	{
		return listOfNodes;
	}
	int getNodeNumber()
	{
		return numberOfNodes;
	}
	int getEdgeNumber()
	{
		int tempCounter=0;
		for(int i=0 ; i<listOfNodes.size() ; i++ )
		{
			tempCounter+=listOfNodes.at(i).getAdjEdgeList().size();
		}
	return tempCounter;
	}
	int getIndexOfid(int idref)
	{
		for(int k=0; k<listOfNodes.size(); k++)
			if(listOfNodes.at(k).getId()==idref)
			{
				return k;
			}
	}

	int isAdjacent(Node a,Node b)
	{
		for(int i=0 ; i<listOfNodes.size() ; i++ )//iterate the nodes
		{
			if( (listOfNodes.at(i).getId()) == a.getId()){//if we find the Node

				for(int j=0; j< listOfNodes.at(j).getAdjEdgeList().size() ; j++)//iterate the list of edges
				{
					if(listOfNodes.at(i).getAdjEdgeList().at(j).getSecondNode()==a.getId())
						return 1;
					else
						return 0;
				}
			}
			else if(listOfNodes.at(i).getId()==b.getId())
			{
				for(int j=0; j< listOfNodes.at(j).getAdjEdgeList().size() ; j++)//iterate the list of edges
				{
					if( listOfNodes.at(i).getAdjEdgeList().at(j).getSecondNode()==b.getId())
							return 1;
						else
							return 0;
				}
			}
		}
		return 0;
	}
	vector<int>& getNeigbors(Node x)
	{

		vector<int> temp;
		for (int i= 0 ; i < listOfNodes.size() ;i++ )
		{
			if(listOfNodes.at(i).getId()==x.getId())//go to node x
			{
				for(int j=0 ;j<listOfNodes.at(i).getAdjEdgeList().size();j++)//iterate the edges
				{
					temp.push_back(listOfNodes.at(i).getAdjEdgeList().at(j).getSecondNode());
				}
			}
		}
		return temp;
	}

	/*
	void add(Node x,Node y,double weight)
	{
		int nodeXIndex=0;
		int nodeYIndex=0;
		int found=0;
		for (int i= 0 ; i < listOfNodes.size() ;i++ )
		{
			if(listOfNodes.at(i).getId()==x.getId())//go to node x
			{
				nodeXIndex=i;
				for(int j=0 ;j<listOfNodes.at(i).getAdjEdgeList().size();j++)//iterate the edges
				{
					if((listOfNodes.at(i).getAdjEdgeList().at(j).getSecondNode()==y.getId()))//check if y is neighbour
					{
						
						if(listOfNodes.at(i).getAdjEdgeList().at(j).getWeight()==weight)
						{
							found=j; //there is an edge so 
							break;
						}
						else
						{
							cout << "Found edge but wrong weight  " << endl;
						}
					}
				}
			}
		}
		if(found ==0)
		{
			for (int i= 0 ; i < listOfNodes.size() ;i++ )
			{
				if(listOfNodes.at(i).getId()==y.getId())//go to node y
				{
					nodeYIndex=i;
				}
			}
			Edge anEdge(y.getId(),x.getId(),weight);
			listOfNodes.at(nodeXIndex).addAdjEdge(anEdge);//add the edges
			Edge anOpEdge(x.getId(),y.getId(),weight);
			listOfNodes.at(nodeYIndex).addAdjEdge(anOpEdge);//to both of the nodes undirected Graph
		}
	}
	*/
	
	double getEdgeValue(Edge& x)
	{
		for(int i=0 ; i <listOfNodes.size() ; i++  )
		{
			for (int j=0 ; j < listOfNodes.at(i).getAdjEdgeList().size() ; j++ )
			{
				if(listOfNodes.at(i).getAdjEdgeList().at(j).getFirstNode()==x.getFirstNode())
					return listOfNodes.at(i).getAdjEdgeList().at(j).getWeight();
			}
		}
	}
	void displayGraph()
	{
		for (int i =0 ; i < listOfNodes.size() ; i++)
		{
			for (int j=0; j < listOfNodes.at(i).getAdjEdgeList().size() ; j++)
			{

				cout << "(" << listOfNodes.at(i).getAdjEdgeList().at(j).getFirstNode() << "->" << listOfNodes.at(i).getAdjEdgeList().at(j).getSecondNode() <<  " cost " << listOfNodes.at(i).getAdjEdgeList().at(j).getWeight()  << ")";
			}
			cout<< " " << endl;
		}
	}
	void displayNodes()
	{
		for (int i =0 ; i < listOfNodes.size() ; i++)
		{
			cout << "Node Number " << listOfNodes.at(i).getId() << endl ;
		}
	}
	void inline generateRandomGraph(double density,double weightFrom,double weightTo)
	{
		//listOfNodes.clear();
		//numberOfNodes=0;
		srand(clock());
		int numberofEdges=(density*(numberOfNodes*(numberOfNodes-1)) )/2;
		while(numberofEdges>0)
		{
			int originating=get_randomNumber(0,numberOfNodes);
			int destination=get_randomNumber(0,numberOfNodes);
			double weight=fRand(weightFrom,weightTo);
			if(originating!=destination){
				Edge edge(originating,destination,weight);
				int temp=0;
				temp=listOfNodes.at(originating).addAdjEdge(edge);
				Edge edgeTwo(destination,originating,weight);
				listOfNodes.at(destination).addAdjEdge(edgeTwo);
				if (temp)
					numberofEdges--;
			}
		}
	}
	inline int get_randomNumber(int rangefrom,int rangeTo)
	{
		return (rand() % rangeTo + rangefrom);
	}
	double fRand(double min, double max)
	{
	    return  (max - min) * ( (double)rand() / (double)RAND_MAX );
	}
	void clearVisited()
	{
		for (int i=0;i<listOfNodes.size();i++)
		{
			for(int j=0; j < listOfNodes.at(i).getAdjEdgeList().size() ;j++)
			{
				listOfNodes.at(i).getAdjEdgeList().at(j).setVisited(false);
			}
		}
	}
	int getNumberOfNodes()
	{
		return numberOfNodes;
	}
};


class ShortestPath{
private:
	Graph* thisGraph;
	std::priority_queue<distanceToNode,std::vector<distanceToNode>,CompdistanceToNode> pqueue;
	vector<distanceToNode> table;
public:
	ShortestPath();
	ShortestPath(Graph* testGraph)
	{
		thisGraph=testGraph;
	}
	void inline DijkstraSSP(Node& startingNode)
		{
			//cout << "inside the shortest path"<< endl;
			int const n=thisGraph->getNumberOfNodes();
			for(int j=0 ; j < thisGraph->getListOfNodes().size() ;j++ )
			{
				distanceToNode anAssoc;
				if(thisGraph->getListOfNodes().at(j).getId()==startingNode.getId()){
					anAssoc.setDistance(0.0);
					anAssoc.setPredecessor(startingNode.getId());
					anAssoc.addCurrent(thisGraph->getListOfNodes().at(j));
					table.push_back(anAssoc);
				}
				else
				{
					anAssoc.setDistance(1000.0);
					anAssoc.addCurrent(thisGraph->getListOfNodes().at(j));
					anAssoc.setPredecessor(NULL);
					table.push_back(anAssoc);
				}
			}
			std::set<Edge> visitedEdges;
			for(int i=0;i<table.size();i++)
				pqueue.push(table.at(i));
			std::set<Node> setOfNodes;
			setOfNodes.insert(startingNode);
			distanceToNode currentLowestDistNode;
			currentLowestDistNode =pqueue.top();
			pqueue.pop();
			while(!pqueue.empty())
			{
				for(int k=0 ; k < currentLowestDistNode.getCurrent().getAdjEdgeList().size() ; k++ )//iterate the edges
				{
					int currentIndex=thisGraph->getIndexOfid(currentLowestDistNode.getCurrent().getId());
					if(!(visitedEdges.find(thisGraph->getListOfNodes()[currentIndex].getAdjEdgeList().at(k))!=visitedEdges.end()))
					{
						double weight;
						weight=currentLowestDistNode.getCurrent().getAdjEdgeList().at(k).getWeight();
						int index;
						index=currentLowestDistNode.getCurrent().getAdjEdgeList().at(k).getSecondNode();//the index of the destination of the Edge
						int currentIndexNode=currentLowestDistNode.getCurrent().getId();
						if( table[currentIndexNode].getDistance()+weight < table.at(index).getDistance() )
						{
							table[index].setDistance( table.at(currentIndexNode).getDistance()+weight);
							int predecessor;
							predecessor=currentLowestDistNode.getCurrent().getId();
							table[index].setPredecessor(predecessor);
							visitedEdges.insert(thisGraph->getListOfNodes()[currentIndex].getAdjEdgeList().at(k));
						}
					}
				}
				currentLowestDistNode=pqueue.top();
				setOfNodes.insert(currentLowestDistNode.getCurrent());
				pqueue.pop();
			}
			thisGraph->clearVisited();
		}
	void printTable()
	{
		for (int tt=0;tt<table.size();tt++)
			cout<< "table " << "Predecessor " << table.at(tt).getPredecessor() << "current" << table.at(tt).getCurrent().getId()<< "table distance " << table.at(tt).getDistance() << endl ;
	}
	vector<Node> path(Node from,Node dest)
	{
		vector<Node> path;
		if (from.getId()!=dest.getId()){//if the source is the same with the destination
			table.clear();
			DijkstraSSP(dest);
			int pred=from.getId();
			double dist;
			while( pred!=dest.getId() && table[pred].getDistance() < 999.9 )
			{
				path.push_back(table[pred].getCurrent());
				pred=table[pred].getPredecessor();
			}
			if(table[pred].getDistance() > 999.9 ){
				cout<< "No Path from "<<from.getId() << " To Node -->" << dest.getId() << endl << endl;;
			}
			else
			{
				path.push_back(table[pred].getCurrent());
				cout << "The path from " << from.getId() << " To Node " << dest.getId() << " is ";
				for(int g=0; g<path.size() ;g++)
				{
					if (g==path.size()-1)
						cout << path.at(g).getId() << endl;
					else
						cout << path.at(g).getId() << "->" ;
				}
			}
			cout << endl;
			table.clear();
			return path;
		}else{
			path.push_back(from);
			table.clear();
			return path;
		}
	}

	inline double path_size(Node from,Node To)
	{
		table.clear();
		DijkstraSSP(To);
		if (table[from.getId()].getDistance() < 999.9)
		{
			 return table[from.getId()].getDistance();
		}
		else
		{
			 return 0.0;
		}
	}

	inline double averageAllPathsToNodeOne(Node destNode)
	{
		double sum=0.0;
		int counter=0;
		DijkstraSSP(destNode);//to create the table
		for(int i=1;i<thisGraph->getListOfNodes().size();i++)
		{
			if(table[i].distance<999.9)
			{
				sum+=table[i].distance;
				counter++;
			}
		}
		return sum/counter;
	}

};

int main() {
	double averageTrials=0.0;
	clock_t startTime = clock();
	for(int i=0;i<10000;i++)//0.20 density
	{
		Graph aGraph(500);
		aGraph.generateRandomGraph(0.20,1.00,10.0);
		ShortestPath aShortestPathClass(&aGraph);
		averageTrials+=aShortestPathClass.averageAllPathsToNodeOne(aGraph.getListOfNodes().at(0));
		//cout << "Average trial number " << i << " an average "<<averageTrials/(i+1) << endl;
	}
	cout <<"The average 1000 trials is: " << averageTrials/10000 <<endl;
	averageTrials=0.0;
	for(int i=0;i<10000;i++)//0.40 density
	{
		Graph aGraph(500);
		aGraph.generateRandomGraph(0.40,1.00,10.0);
		ShortestPath aShortestPathClass(&aGraph);
		averageTrials+=aShortestPathClass.averageAllPathsToNodeOne(aGraph.getListOfNodes().at(0));
		//cout << "Average trial number " << i << " an average "<<averageTrials/i << endl;
	}
	cout <<"Graph density 0.40 The average 1000 trials is: " << averageTrials/10000 <<endl;
	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
	return 0;
}
