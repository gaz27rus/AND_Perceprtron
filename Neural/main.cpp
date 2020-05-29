#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

class neuron
{
private:
	double input, output, delta, GRAD; 
public:
	double weight_delta = 0;
	double weight;
	void set_weight(double Weight);
	void set_input(double num);
	void set_output(double num);
	void set_delta(double num);
	void set_GRAD(double num);
	double get_input();
	double get_GRAD();
	double get_output();
	double get_weight();
	double normalized();
	double brain_normalized(int index);
	double get_delta();
	vector<double> Weight;
	vector<double> grad;
	vector<double> Weight_Delta;
};

class Neural_Network
{
public:
	Neural_Network(int size_of_network, double Speed, double Moment, int hidden_count) : learningSpeed(Speed), moment(Moment)
	{
		brain.resize(size_of_network);

		for (int i = 0; i < size_of_network; i++)
		{
			brain[i].Weight.resize(hidden_count);
			brain[i].grad.resize(hidden_count);
			brain[i].Weight_Delta.resize(hidden_count);
			for (int j = 0; j < hidden_count; j++)
			{
				brain[i].Weight[j] = (DoubleRand(1.0, 0));
			}
			
			//cout << brain[i].get_weight() << endl;
		}

		hidden.resize(hidden_count);

		for (int i = 0; i < size_of_network; i++)
		{
			hidden[i].set_weight(DoubleRand(1.0, 0));
			//cout << brain[i].get_weight() << endl;
		}

		//hidden.set_weight(DoubleRand(1.0, 0));
		input_bias.set_input(1.0);
		input_bias.set_output(1.0);
		input_bias.Weight.resize(hidden_count);
		input_bias.grad.resize(hidden_count);
		input_bias.Weight_Delta.resize(hidden_count);
		for (int i = 0; i < hidden_count; i++)
		{
			input_bias.Weight[i] = (DoubleRand(1.0, 0));
		}
		
		//cout << hidden.get_weight() << endl;
		//cout << input_bias.get_weight() << endl;
	};

	
	double DoubleRand(double _max, double _min);
	neuron* get_pointer(int index);
	neuron* get_hidden_pointer(int index);
	neuron* get_bias_pointer();
	neuron* get_output_pointer();
	double normalized_sum(int ind);
	double normalized_hidden_sum();
	void count_error(double num, double expect);
	double get_error();
	void pass_error();
	void activation(neuron* hidden, double num);
	void output_activation(neuron* output, double num);
	void Backpropagation(double out_ideal, double out_actual);

private:
	//neuron hidden;
	neuron input_bias;
	neuron output;
	vector<neuron> brain;
	vector<neuron> hidden;
	//vector<double> it_error;
	double moment;
	double learningSpeed;
	double error;
};



int main()
{
	srand(time(NULL));
	int size, hidden_size;
	double speed, moment;
	cout << "Enter count of input neurons: "; cin >> size;
	cout << "Enter count of hidden neurons: "; cin >> hidden_size;
	cout << "Enter speed of learning: "; cin >> speed;
	cout << "Enter moment of learning: "; cin >> moment;
	Neural_Network  perceptron(size, speed, moment, hidden_size);

	for (int i = 0; i < 10000; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//cout << "Set " << j + 1 << endl;
			/*if (j == 0)
			{
				for (int i = 0; i < size; i++)
				{
					int a;
					//cin >> a;
					perceptron.get_pointer(i)->set_input(1);
					perceptron.get_pointer(i)->set_output(1);
				}

				perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout << perceptron.get_output_pointer()->get_output() << endl;


				int k = 0;
				/*cout << "Enter the answer: ";
				cin >> k;

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
				/*cout << "New weights: " << endl;

				for (int i = 0; i < size; i++)
				{
					cout << perceptron.get_pointer(i)->get_weight() << endl;
				}
				cout << perceptron.get_hidden_pointer()->get_weight() << endl;
				cout << perceptron.get_bias_pointer()->get_weight() << endl;
			}
			if (j == 1)
			{
				for (int i = 0; i < size; i++)
				{
					int a;
					//cin >> a;
					perceptron.get_pointer(i)->set_input(0);
					perceptron.get_pointer(i)->set_output(0);
				}

				perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout << perceptron.get_output_pointer()->get_output() << endl;


				int k = 0;
				/*cout << "Enter the answer: ";
				cin >> k;

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
				/*cout << "New weights: " << endl;

				for (int i = 0; i < size; i++)
				{
				cout << perceptron.get_pointer(i)->get_weight() << endl;
				}
				cout << perceptron.get_hidden_pointer()->get_weight() << endl;
				cout << perceptron.get_bias_pointer()->get_weight() << endl;
			}*/
			if (j == 0)
			{
				perceptron.get_pointer(0)->set_input(1);
				perceptron.get_pointer(0)->set_output(1);
				perceptron.get_pointer(1)->set_input(0);
				perceptron.get_pointer(1)->set_output(0);

				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.get_hidden_pointer(k)->set_input(perceptron.normalized_sum(k));
				}
				//perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.activation(perceptron.get_hidden_pointer(k), perceptron.get_hidden_pointer(k)->get_input());
				}
				//perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_output_pointer()->get_input());
				//perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout <<"1 0: " << perceptron.get_output_pointer()->get_output() << endl;


				int k = 0;
				/*cout << "Enter the answer: ";
				cin >> k;*/

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
				/*cout << "New weights: " << endl;

				for (int i = 0; i < size; i++)
				{
				cout << perceptron.get_pointer(i)->get_weight() << endl;
				}
				cout << perceptron.get_hidden_pointer()->get_weight() << endl;
				cout << perceptron.get_bias_pointer()->get_weight() << endl;*/
			}
			if (j == 1)
			{
				perceptron.get_pointer(0)->set_input(0);
				perceptron.get_pointer(0)->set_output(0);
				perceptron.get_pointer(1)->set_input(1);
				perceptron.get_pointer(1)->set_output(1);

				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.get_hidden_pointer(k)->set_input(perceptron.normalized_sum(k));
				}
				//perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.activation(perceptron.get_hidden_pointer(k), perceptron.get_hidden_pointer(k)->get_input());
				}
				//perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_output_pointer()->get_input());
				//perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout << "0 1: " << perceptron.get_output_pointer()->get_output() << endl;


				int k = 0;
				/*cout << "Enter the answer: ";
				cin >> k;*/

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
			}
			if (j == 2)
			{
				perceptron.get_pointer(0)->set_input(1);
				perceptron.get_pointer(0)->set_output(1);
				perceptron.get_pointer(1)->set_input(1);
				perceptron.get_pointer(1)->set_output(1);

				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.get_hidden_pointer(k)->set_input(perceptron.normalized_sum(k));
				}
				//perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.activation(perceptron.get_hidden_pointer(k), perceptron.get_hidden_pointer(k)->get_input());
				}
				//perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_output_pointer()->get_input());
				//perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout << "1 1: " << perceptron.get_output_pointer()->get_output() << endl;


				int k = 1;
				/*cout << "Enter the answer: ";
				cin >> k;*/

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
				/*cout << "New weights: " << endl;

				for (int i = 0; i < size; i++)
				{
				cout << perceptron.get_pointer(i)->get_weight() << endl;
				}
				cout << perceptron.get_hidden_pointer()->get_weight() << endl;
				cout << perceptron.get_bias_pointer()->get_weight() << endl;*/
			}
			if (j == 3)
			{
				perceptron.get_pointer(0)->set_input(0);
				perceptron.get_pointer(0)->set_output(0);
				perceptron.get_pointer(1)->set_input(0);
				perceptron.get_pointer(1)->set_output(0);

				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.get_hidden_pointer(k)->set_input(perceptron.normalized_sum(k));
				}
				//perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
				for (int k = 0; k < hidden_size; k++)
				{
					perceptron.activation(perceptron.get_hidden_pointer(k), perceptron.get_hidden_pointer(k)->get_input());
				}
				//perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
				perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
				perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_output_pointer()->get_input());
				//perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
				//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
				cout << "0 0: " << perceptron.get_output_pointer()->get_output() << endl;


				int k = 0;
				/*cout << "Enter the answer: ";
				cin >> k;*/

				perceptron.count_error(k, perceptron.get_output_pointer()->get_output());
				perceptron.Backpropagation(k, perceptron.get_output_pointer()->get_output());
				/*cout << "New weights: " << endl;

				for (int i = 0; i < size; i++)
				{
				cout << perceptron.get_pointer(i)->get_weight() << endl;
				}
				cout << perceptron.get_hidden_pointer()->get_weight() << endl;
				cout << perceptron.get_bias_pointer()->get_weight() << endl;*/
			}
		}
		//perceptron.Backpropagation()
		cout << endl << "Error: " << ((perceptron.get_error() / 4) * 100) << "%"  << endl << endl;
		perceptron.pass_error();
	}
	/*for (int i = 0; i < size; i++)
	{
		int a;
		cin >> a;
		perceptron.get_pointer(i)->set_input(a);
		perceptron.get_pointer(i)->set_output(a);
	}

	perceptron.get_hidden_pointer()->set_input(perceptron.normalized_sum());
	perceptron.activation(perceptron.get_hidden_pointer(), perceptron.get_hidden_pointer()->get_input());
	perceptron.get_output_pointer()->set_input(perceptron.normalized_hidden_sum());
	perceptron.output_activation(perceptron.get_output_pointer(), perceptron.get_hidden_pointer()->get_output());
	//perceptron.get_output_pointer()->set_output(perceptron.get_hidden_pointer()->get_output());
	cout << perceptron.get_output_pointer()->get_output() << endl;*/

	system("pause");
	return(0);
}


double Neural_Network::DoubleRand(double _max, double _min)
{
	return _min + double(rand()) / RAND_MAX * (_max - _min);
}

void neuron::set_input(double num)
{
	input = num;
}

void neuron::set_weight(double Weight)
{
	weight = Weight;
}

neuron* Neural_Network::get_pointer(int index)
{
	return &(brain[index]);
}

double neuron::normalized()
{
	return input * weight;
}

double neuron::brain_normalized(int index)
{
	return input * Weight[index];
}

neuron* Neural_Network::get_hidden_pointer(int index)
{
	return &(hidden[index]);
}

neuron* Neural_Network::get_bias_pointer()
{
	return &(input_bias);
}

double Neural_Network::normalized_hidden_sum()
{
	double sum = 0;
	for (int l = 0; l < hidden.size(); l++)
	{
		sum += hidden[l].normalized();
	}
	return sum;
}

double Neural_Network::normalized_sum(int ind)
{
	double sum = 0;
	for (int i = 0; i < brain.size(); i++)
	{
		sum += brain[i].brain_normalized(ind);
	}
	sum += input_bias.brain_normalized(ind);
	return sum;
}

void Neural_Network::activation(neuron* hidden, double num)
{
	(*hidden).set_output((1 / (1 + exp(-num))));
}

void Neural_Network::output_activation(neuron* output, double num)
{
	(*output).set_output((1 / (1 + exp(-num))));
}

void neuron::set_output(double num)
{
	output = num;
}

double neuron::get_input()
{
	return input;
}

void neuron::set_delta(double num)
{
	delta = num;
}

void neuron::set_GRAD(double num)
{
	GRAD = num;
}

neuron* Neural_Network::get_output_pointer()
{
	return &(output);
}

double neuron::get_output()
{
	return output;
}

double neuron::get_GRAD()
{
	return GRAD;
}

double neuron::get_delta()
{
	return delta;
}

double neuron::get_weight()
{
	return weight;
}

void Neural_Network::count_error(double num, double expect)
{
	error += pow((num - expect), 2.0);// / brain.size();
	//cout << "Error: " << error << endl;
}

double Neural_Network::get_error()
{
	return error;
}

void Neural_Network::pass_error()
{
	error = 0;
}

void Neural_Network::Backpropagation(double out_ideal, double out_actual)
{
	output.set_delta((out_ideal - out_actual) * ((1 - out_actual)*out_actual)); // дельта выходного

	for (int l = 0; l < hidden.size(); l++)
	{
		hidden[l].set_delta(((1 - hidden[l].get_output()) * hidden[l].get_output()) * (hidden[l].get_weight() * output.get_delta())); //дельта скрытого
		hidden[l].set_GRAD(hidden[l].get_output() * output.get_delta());
		hidden[l].set_weight(hidden[l].get_weight() + (learningSpeed * hidden[l].get_GRAD() + moment * hidden[l].weight_delta)); //поменять с учетом изменения веса и момента
		hidden[l].weight_delta = learningSpeed * hidden[l].get_GRAD() + moment * hidden[l].weight_delta;
	}

	for (int i = 0; i < brain.size(); i++)
	{
		//cout << brain[i].get_GRAD();
		for (int l = 0; l < brain[i].Weight.size(); l++)
		{
			brain[i].grad[l] = brain[i].get_output() * hidden[l].get_delta();
			brain[i].Weight[l] = brain[i].Weight[l] + (learningSpeed * brain[i].grad[l] + moment * brain[i].Weight_Delta[l]);
			brain[i].Weight_Delta[l] = learningSpeed * brain[i].grad[l] + moment * brain[i].Weight_Delta[l];
			//brain[i].set_weight(brain[i].get_weight() + (learningSpeed * brain[i].get_GRAD() + moment * brain[i].weight_delta));
			//brain[i].set_GRAD(brain[i].get_output() * hidden.get_delta());
		}
		//brain[i].set_GRAD(brain[i].get_output() * hidden.get_delta());
		//brain[i].set_weight(brain[i].get_weight() + (learningSpeed * brain[i].get_GRAD() + moment * brain[i].weight_delta)); //поменять с учетом изменения веса и момента
		//brain[i].weight_delta = learningSpeed * brain[i].get_GRAD() + moment * brain[i].weight_delta;
	}
	for (int l = 0; l < input_bias.Weight.size(); l++)
	{
		input_bias.grad[l] = input_bias.get_output() * hidden[l].get_delta();
		input_bias.Weight[l] = input_bias.Weight[l] + (learningSpeed * input_bias.grad[l] + moment * input_bias.Weight_Delta[l]);
		input_bias.Weight_Delta[l] = learningSpeed * input_bias.grad[l] + moment * input_bias.Weight_Delta[l];
	}
	/*input_bias.set_GRAD(input_bias.get_output() * hidden.get_delta());
	input_bias.set_weight(input_bias.get_weight() + (learningSpeed * input_bias.get_GRAD() + moment * input_bias.weight_delta));
	input_bias.weight_delta = learningSpeed * input_bias.get_GRAD() + moment * input_bias.weight_delta;*/
}