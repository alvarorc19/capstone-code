#ifndef MODELBASE_H
#define MODELBASE_H

class ModelBase {
    public:
        virtual ~ModelBase() = default;
        virtual void change_spin_randomly(int index) = 0;
        virtual void change_spin_randomly(ivec indices) = 0;
        virtual double compute_magnetisation() = 0;
        virtual double compute_single_partition_function() = 0;
        virtual void write_magnetisation() = 0;
        virtual void write_energy() = 0;
        
};

#endif
