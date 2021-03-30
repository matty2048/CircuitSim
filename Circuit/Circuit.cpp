// Circuit.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <Eigen/Dense>
#include <eigen3/unsupported/Eigen/FFT>
#include <math.h>
#include <complex>


#include <implot.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <gl/glew.h>
#include <GLFW/glfw3.h>


using namespace Eigen;
typedef Eigen::dcomplex jfloat;

struct Node
{
	Node():ID(0),V(0) {};
	Node(char* ID) : ID(ID) {};
	const char* ID;
	double V;
};
enum class Type
{
	Resistor,
	Capacitor,
	Inductor,
	Vsource,
	Isource
};
class Component
{
public:
	Component() {};
	//Component(int ID, double R) : ID(ID), G(1/R) {};
	int ID;
	int nodeA;
	int nodeB;
	jfloat G; //conductivity
	jfloat GetG(double freq) {
		return this->GetConductance(freq);
	}
	virtual jfloat GetConductance(double freq) { return 0; };
	Type type;
};

class Resistor :public Component
{
public:
	Resistor(int ID, double R) {
		this->ID = ID;
		this->G = 1 / R;
		this->type = Type::Resistor;
	};
	virtual jfloat GetConductance(double freq) override
	{
		return G;
	} 
};

class Capacitor : public Component
{
public:
	Capacitor(int ID, double C) {
		this->ID = ID;
		this->C = C;
		this->type = Type::Capacitor;
	}
	double C;
	virtual jfloat GetConductance(double freq) override
	{
		if (freq == 0) return 0;
		return std::complex<double>(0,1) * 2.0 * double(EIGEN_PI) * double(freq) * C;
	}
};

class Inductor : public Component
{
public:
	Inductor(int ID, double L) {
		this->ID = ID;
		this->L = L;
		this->type = Type::Inductor;
	};
	double L;
	
	jfloat GetConductance(double freq) override
	{
		if (freq == 0) return 1;
		return (std::complex<double>(1,0) / (std::complex<double>(0,1) * 2.0 * double(EIGEN_PI) * freq * L));
	}
};

struct VoltageSource : public Component
{

	double V;
	int Link = 0;
	VoltageSource(int ID, double V)  {
		this->ID = ID;
		this->V = V;
		this->type = Type::Vsource;
	};
};

struct CurrentSource : public Component
{
	double I;
	int Link = 0; //if this source is linked to a Geq
	CurrentSource(int ID, double I)
	{
		this->ID = ID;
		this->I = I;
		this->type = Type::Isource;
	};

};


class Graph
{
public:
	Graph() {};
	int n_sources = 0;
	std::vector<Node*> Nodes;
	std::vector<Component*> Connections;
	std::vector<VoltageSource*> sources;
	std::vector<CurrentSource*> Csources;
	unsigned int operator [](const char* ID) //returns a node from an NodeID
	{
		unsigned int i = 0;
		Node* n = Nodes[i];
		if (n->ID == ID) return 0;
		while (n->ID != ID && i < Nodes.size())
		{
			n = Nodes[i];
			i++;
		}
		return i-1;
	}
	unsigned int operator [](const int ID) //returns a connection from connection ID
	{
		unsigned int i = 0;
		Component* c = Connections[i];
		if (c->ID == ID) return 0;
		while (c->ID != ID && i < Connections.size())
		{
			c = Connections[i];
			i++;
		}
		return i-1;
	}
	void AddNode(const char* ID)
	{
		Node* n = new Node();
		n->ID = ID;
		Nodes.emplace_back(n);
	}
	
	void AddResistor(int ID, double R, const char* NodeA, const char* NodeB)
	{
		Component* c = new Resistor(ID, R);
		Connections.push_back(c);
		Connections.back()->nodeA = this->operator[](NodeA);
		Connections.back()->nodeB = this->operator[](NodeB);
	}
	void AddCapacitor(int ID, double C, const char* NodeA, const char* NodeB)
	{
		Component* c = new Capacitor(ID, C);
		Connections.push_back(c);
		Connections.back()->nodeA = this->operator[](NodeA);
		Connections.back()->nodeB = this->operator[](NodeB);
	}
	void AddInductor(int ID, double L, const char* NodeA, const char* NodeB)
	{
		Component* c = new Inductor(ID, L);
		Connections.push_back(c);
		Connections.back()->nodeA = this->operator[](NodeA);
		Connections.back()->nodeB = this->operator[](NodeB);
	}

	void AddSource(int ID, double I, const char* NodeA, const char* NodeB)
	{
		VoltageSource* c = new VoltageSource(ID,I);
		sources.push_back(c);
		sources.back()->nodeA = this->operator[](NodeA);
		sources.back()->nodeB = this->operator[](NodeB);
	}

	void AddCsource(int ID, double I, const char* NodeA, const char* NodeB)
	{
		CurrentSource* c = new CurrentSource(ID, I);
		Csources.push_back(c);
		Csources.back()->nodeA = this->operator[](NodeA);
		Csources.back()->nodeB = this->operator[](NodeB);
	}

	int n() {return Nodes.size() - 1; };
	int m() { return sources.size(); };
	jfloat Get_Sum(int a, double freq) //gets the sum of conductances connected to a node
	{
		jfloat res = 0;
		for (unsigned int i = 0; i < Connections.size(); i++)
		{
			if (Connections[i]->nodeA == a | Connections[i]->nodeB == a )
			{
				res = res + Connections[i]->GetG(freq);
			}
		}
		return res;
	}
	jfloat FindConnections(int a, int b, double freq)
	{
		
		for (unsigned int i = 0; i < Connections.size(); i++)
		{
			if ((((Connections[i]->nodeA == a) && (Connections[i]->nodeB == b)) || ((Connections[i]->nodeA == b) && (Connections[i]->nodeB == a))))
			{
				return -Connections[i]->GetG(freq);
			}
		}
		return 0;

	}
	void PrintMatrix(MatrixXcd matrix)
	{
		std::cout << std::endl;
		for (unsigned int i = 0; i < matrix.rows(); i++)
		{
			for (unsigned int j = 0; j < matrix.cols(); j++)
				std::cout << matrix(i,j);
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	void Calculate(double freq)
	{
		MatrixXcd A_matrix(n() + m(), n() + m());
		MatrixXcd G_matrix(n(), n());
		MatrixXcd B_matrix(n(), m());
		MatrixXcd C_matrix(m(), n());
		MatrixXcd D_matrix(m(), m());


		//creates G matrix
		for (unsigned int i = 0; i < n(); i++)
		{
			for (unsigned int j = 0; j < n(); j++)
			{
				if (j == i) G_matrix(j, i) = Get_Sum(j+1,freq);
				else G_matrix(j, i) = FindConnections(i+1, j+1,freq);
			}
		}
		//PrintMatrix(G_matrix);
		//creates B matrix
		for (unsigned int i = 0; i < m(); i++) //loops over every voltage source
		{
			for (unsigned int j = 1; j < n() + 1; j++) //loops over every node
			{
				B_matrix(j - 1, i) = 0;
				if (this->operator[](Nodes[j]->ID) == sources[i]->nodeA) //if + terminal is connected to the node 
				{
					B_matrix(j - 1, i) = 1;
				}
				if (this->operator[](Nodes[j]->ID) == sources[i]->nodeB) //if + terminal is connected to the node 
				{
					B_matrix(j - 1, i) = -1;
				}
			}
		}
		
		//PrintMatrix(B_matrix);
		
		//creates C matrix
		C_matrix = B_matrix.transpose();
		D_matrix.setZero();

		A_matrix << G_matrix , B_matrix, C_matrix, D_matrix;
		//for (int i = 0; i < n() + m(); i++)
		//{
		//	for (unsigned int j = 0; j < n() + m(); j++)
		//		std::cout << A_matrix(i, j) << " ";
		//	std::cout << std::endl;
		//}

		MatrixXd Z_matrix(m() + n(), 1);
		Z_matrix.setZero();
		
		for (unsigned int i = 0; i< n(); i++) //loops over all nodes
		{
			for (auto s : Csources)
			{
				if (s->nodeA == i + 1)
				{
					Z_matrix(i) += s->I;
				}
				if (s->nodeB == i + 1)
				{
					Z_matrix(i) -= s->I;
				}
			}
		}
		for (unsigned int i = 0; i < m(); i++)
		{
			Z_matrix((m() + n()) - i - 1 , 0) = sources[i]->V;
		}


		X_matrix = A_matrix.inverse() * (Z_matrix);

		for (unsigned int j = 0; j < n() + m(); j++) {
			if (j < n())std::cout << "Node:" << Nodes[j + 1]->ID << " ";
			else std::cout <<"Current from:" << sources[j - n()]->ID << " ";
			std::cout << abs(X_matrix(j, 0));
			if (j < n())std::cout << "V";
			else std::cout << "A";
			std::cout << std::endl;
		}
	}

	MatrixXcd X_matrix;
	void TranSolve(float sT, float eT, float dT)
	{
		for (unsigned int i = 0; i < Connections.size(); i++)
		{
			if (Connections[i]->type == Type::Capacitor) //replaces the capacitor with a current source and resistor
			{
				Capacitor* c = dynamic_cast<Capacitor*>(Connections[i]);
				Connections.erase(Connections.begin() + i); //removes the component

				this->AddResistor(Connections.size(), c->C / dT, Nodes[c->nodeA]->ID, Nodes[c->nodeB]->ID);
				this->AddCsource(Csources.size(), 0, Nodes[c->nodeA]->ID, Nodes[c->nodeB]->ID);
				Csources.back()->Link = this->operator[](Connections.back()->ID); //links the current and resistor
			}
			if (Connections[i]->type == Type::Inductor)
			{
				Inductor* c = dynamic_cast<Inductor*>(Connections[i]);
				Connections.erase(Connections.begin() + i);

				this->AddResistor(Connections.size(), c->L / dT, Nodes[c->nodeA]->ID, Nodes[c->nodeB]->ID);
				this->AddCsource(Connections.size(), 0, Nodes[c->nodeB]->ID, Nodes[c->nodeA]->ID);
			}

		}
		X_matrix = MatrixXcd(m() + n(), 1).setZero();
		
	}
	void TranStep()
	{
		this->Calculate(0);
		for (unsigned int i = 0; i < Csources.size(); i++)
		{
			if (Csources[i]->Link != 0)
			{
				std::complex<double> T = X_matrix(Csources[i]->nodeA - 1);
				Csources[i]->I = (Connections[Csources[i]->Link]->GetConductance(0) * T).real();
			}
			PrintMatrix(X_matrix);
		}
	}
};



int main()
{


	if (!glfwInit())
		return -1;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);  //initializes GLFW
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(1000, 1000, "ViewPort", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}
	//Creates GLFW window for drawing to
	glfwMakeContextCurrent(window);

	if (glewInit() != GLEW_OK) ///initializes glew
	{
		std::cout << "error with glew :c" << std::endl;
		std::cin;
		return -1;
	}

	glfwSwapInterval(1);


	glEnable(GL_DEBUG_OUTPUT);
	glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);       //enables glError + binds message callback on error




	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	ImGuiIO& io = ImGui::GetIO();
	// Setup Platform/Renderer bindings
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 330");
	// Setup Dear ImGui style
	ImGui::StyleColorsDark();



	Graph g;
	g.AddNode("end"); //this is 0v, 0v always first

	g.AddNode("v1");
	g.AddNode("v2");


	g.AddSource(1, 1, "v1", "end");
	g.AddResistor(2, 200, "v1", "v2");
	g.AddCapacitor(3, 0.0001, "v2", "end");





	g.TranSolve(0, 1, 0.0001);

	float T = 0;
	std::vector<double> Tx;
	std::vector<double> Vx;
	Tx.emplace_back(T);
	while (!glfwWindowShouldClose(window))
	{


		glfwPollEvents();
		glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
		glClear(GL_COLOR_BUFFER_BIT);

		// feed inputs to dear imgui, start new frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// rendering our geometries
		g.TranStep();
		Vx.emplace_back(g.X_matrix(1).real());
		T += 0.0001;
		Tx.emplace_back(T);
		if (ImPlot::BeginPlot("graph")) {
			ImPlot::PlotLine("v2", Tx.data(), Vx.data(), Tx.size());
			ImPlot::EndPlot();
		}
		//// render your GUI
		ImGui::Begin("Demo window");
		ImGui::Button("Hello!");
		ImGui::End();
		

		// Render dear imgui into screen
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glfwSwapBuffers(window);
}
	std::cout << "got to the end" << std::endl;
}