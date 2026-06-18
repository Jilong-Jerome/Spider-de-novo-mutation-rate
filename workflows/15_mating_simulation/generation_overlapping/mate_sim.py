import random
import matplotlib.pyplot as plt

# -------------------------------
# Simulation Parameters
# -------------------------------
# Initial numbers in the first generation:
initial_num_males = 50
initial_num_females = 50

# Simulation time parameters:
T_max = 500.0   # total simulation time
dt = 0.1        # time step

# Encounter and death rate parameters:
alpha = 0.0001           # base encounter rate constant (per female per time unit)
death_rate_moving = 0.1   # death rate for a male while moving (per time unit)
death_rate_mating = 0.6    # probability a male dies during a mating event
death_rate_female = 0.01  # baseline death rate for females (per time unit)

# Generation overlapping parameters:
generation_interval = 2.0       # every 50 time units, introduce a new generation
new_generation_males = 50          # number of new males each generation
new_generation_females = 50        # number of new females each generation

# -------------------------------
# Helper functions to create new individuals.
# -------------------------------
males = []
females = []
next_male_id = 0
next_female_id = 0

def create_males(num, current_time):
    """Create a list of new male individuals."""
    global next_male_id
    new_males = []
    for _ in range(num):
        new_males.append({
            'id': next_male_id,
            'birth_time': current_time,         # record when this individual is introduced
            'alive': True,
            'first_mating_time': None,          # will be set when the individual mates for the first time (relative to its birth_time)
            'sex': 'male'
        })
        next_male_id += 1
    return new_males

def create_females(num, current_time):
    """Create a list of new female individuals."""
    global next_female_id
    new_females = []
    for _ in range(num):
        new_females.append({
            'id': next_female_id,
            'birth_time': current_time,         # record when this individual is introduced
            'alive': True,
            'first_mating_time': None,          # will be set when the individual mates for the first time (relative to its birth_time)
            'sex': 'female'
        })
        next_female_id += 1
    return new_females

# -------------------------------
# Initialize the first generation at time 0.
# -------------------------------
males.extend(create_males(initial_num_males, 0))
females.extend(create_females(initial_num_females, 0))

# -------------------------------
# Run the Simulation with Overlapping Generations
# -------------------------------
time = 0.0
next_generation_time = generation_interval  # when to add the next generation

while time < T_max:
    # Introduce a new generation at the specified interval.
    if time >= next_generation_time:
        males.extend(create_males(new_generation_males, time))
        females.extend(create_females(new_generation_females, time))
        next_generation_time += generation_interval

    # Update females (stationary, low death rate).
    for female in females:
        if female['alive']:
            if random.random() < death_rate_female * dt:
                female['alive'] = False

    # Update males (active searchers).
    for male in males:
        if not male['alive']:
            continue  # skip dead males

        # Males may die while moving.
        if random.random() < death_rate_moving * dt:
            male['alive'] = False
            continue

        # Identify all currently alive females.
        living_females = [f for f in females if f['alive']]
        if living_females:
            # Encounter probability scales with the number of living females.
            encounter_probability = alpha * len(living_females) * dt
            if random.random() < encounter_probability:
                # Record the male's first mating time if not already recorded.
                # Record the time elapsed since the male was introduced.
                if male['first_mating_time'] is None:
                    male['first_mating_time'] = time - male['birth_time']
                
                # Randomly choose a female.
                chosen_female = random.choice(living_females)
                # Record the female's first mating time (time elapsed since her introduction) if not already recorded.
                if chosen_female['first_mating_time'] is None:
                    chosen_female['first_mating_time'] = time - chosen_female['birth_time']

                # The cost of mating: the male may die immediately.
                if random.random() < death_rate_mating:
                    male['alive'] = False

    time += dt

# -------------------------------
# Summarize the Results and Write Report Files
# -------------------------------
# Collect first mating times for all individuals that mated.
male_mating_times = [male['first_mating_time'] for male in males if male['first_mating_time'] is not None]
female_mating_times = [female['first_mating_time'] for female in females if female['first_mating_time'] is not None]

print("Male first mating times (relative to birth):")
print(male_mating_times)
print("Female first mating times (relative to birth):")
print(female_mating_times)

# Write simulation parameters to a text file.
with open("simulation_parameters.txt", "w") as f:
    f.write("Simulation Parameters:\n")
    f.write(f"initial_num_males: {initial_num_males}\n")
    f.write(f"initial_num_females: {initial_num_females}\n")
    f.write(f"T_max: {T_max}\n")
    f.write(f"dt: {dt}\n")
    f.write(f"alpha: {alpha}\n")
    f.write(f"death_rate_moving: {death_rate_moving}\n")
    f.write(f"death_rate_mating: {death_rate_mating}\n")
    f.write(f"death_rate_female: {death_rate_female}\n")
    f.write(f"generation_interval: {generation_interval}\n")
    f.write(f"new_generation_males: {new_generation_males}\n")
    f.write(f"new_generation_females: {new_generation_females}\n")

# Write a tab-delimited file with the first mating record for each individual that mated.
with open("first_mating_record.txt", "w") as f:
    f.write("ind_id\tbirth_time\trelative_first_mating_time\tsex\n")
    for male in males:
        if male['first_mating_time'] is not None:
            f.write(f"male_{male['id']}\t{male['birth_time']}\t{male['first_mating_time']}\t{male['sex']}\n")
    for female in females:
        if female['first_mating_time'] is not None:
            f.write(f"female_{female['id']}\t{female['birth_time']}\t{female['first_mating_time']}\t{female['sex']}\n")

# -------------------------------
# Plot the Distribution of Relative First Mating Times
# -------------------------------
plt.figure(figsize=(10, 5))
plt.hist(male_mating_times, bins=30, alpha=0.5, label="Males")
plt.hist(female_mating_times, bins=30, alpha=0.5, label="Females")
plt.xlabel("Relative First Mating Time (time since introduction)")
plt.ylabel("Count")
plt.title("Distribution of Relative First Mating Times (Overlapping Generations)")
plt.legend()
plt.show()

