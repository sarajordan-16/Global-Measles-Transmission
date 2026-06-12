#!/usr/bin/env python3
"""
Replace country names in FASTA headers with their corresponding continents.
Header format: >accession/date/country/length/genotype
Will become: >accession/date/continent/length/genotype
"""

from Bio import SeqIO

# Dictionary mapping countries to continents
COUNTRY_TO_CONTINENT = {
    # Africa
    'Algeria': 'Africa',
    'Angola': 'Africa',
    'Benin': 'Africa',
    'Botswana': 'Africa',
    'Burkina_Faso': 'Africa',
    'Burundi': 'Africa',
    'Cameroon': 'Africa',
    'Cape_Verde': 'Africa',
    'Central_African_Republic': 'Africa',
    'Chad': 'Africa',
    'Comoros': 'Africa',
    'Congo': 'Africa',
    'Democratic_Republic_of_the_Congo': 'Africa',
    "Cote_d'Ivoire": 'Africa',
    'Djibouti': 'Africa',
    'Egypt': 'Africa',
    'Equatorial_Guinea': 'Africa',
    'Eritrea': 'Africa',
    'Ethiopia': 'Africa',
    'Gabon': 'Africa',
    'Gambia': 'Africa',
    'Ghana': 'Africa',
    'Guinea': 'Africa',
    'Guinea-Bissau': 'Africa',
    'Kenya': 'Africa',
    'Lesotho': 'Africa',
    'Liberia': 'Africa',
    'Libya': 'Africa',
    'Madagascar': 'Africa',
    'Malawi': 'Africa',
    'Mali': 'Africa',
    'Mauritania': 'Africa',
    'Mauritius': 'Africa',
    'Morocco': 'Africa',
    'Mozambique': 'Africa',
    'Namibia': 'Africa',
    'Niger': 'Africa',
    'Nigeria': 'Africa',
    'Rwanda': 'Africa',
    'Sao_Tome_and_Principe': 'Africa',
    'Senegal': 'Africa',
    'Seychelles': 'Africa',
    'Sierra_Leone': 'Africa',
    'Somalia': 'Africa',
    'South_Africa': 'Africa',
    'South_Sudan': 'Africa',
    'Sudan': 'Africa',
    'Swaziland': 'Africa',
    'Tanzania': 'Africa',
    'Togo': 'Africa',
    'Tunisia': 'Africa',
    'Uganda': 'Africa',
    'Zambia': 'Africa',
    'Zimbabwe': 'Africa',
    
    # Asia
    'Afghanistan': 'Asia',
    'Armenia': 'Asia',
    'Azerbaijan': 'Asia',
    'Bahrain': 'Asia',
    'Bangladesh': 'Asia',
    'Bhutan': 'Asia',
    'Brunei': 'Asia',
    'Cambodia': 'Asia',
    'China': 'Asia',
    'Georgia': 'Asia',
    'India': 'Asia',
    'Indonesia': 'Asia',
    'Iran': 'Asia',
    'Iraq': 'Asia',
    'Israel': 'Asia',
    'Japan': 'Asia',
    'Jordan': 'Asia',
    'Kazakhstan': 'Asia',
    'Kuwait': 'Asia',
    'Kyrgyzstan': 'Asia',
    'Laos': 'Asia',
    'Lebanon': 'Asia',
    'Malaysia': 'Asia',
    'Maldives': 'Asia',
    'Mongolia': 'Asia',
    'Myanmar': 'Asia',
    'Nepal': 'Asia',
    'North_Korea': 'Asia',
    'Oman': 'Asia',
    'Pakistan': 'Asia',
    'Palestine': 'Asia',
    'Philippines': 'Asia',
    'Qatar': 'Asia',
    'Saudi_Arabia': 'Asia',
    'Singapore': 'Asia',
    'South_Korea': 'Asia',
    'Sri_Lanka': 'Asia',
    'Syria': 'Asia',
    'Taiwan': 'Asia',
    'Tajikistan': 'Asia',
    'Thailand': 'Asia',
    'Timor-Leste': 'Asia',
    'Turkey': 'Asia',
    'Turkmenistan': 'Asia',
    'United_Arab_Emirates': 'Asia',
    'Uzbekistan': 'Asia',
    'Vietnam': 'Asia',
    'Yemen': 'Asia',
    
    # Europe
    'Albania': 'Europe',
    'Andorra': 'Europe',
    'Austria': 'Europe',
    'Belarus': 'Europe',
    'Belgium': 'Europe',
    'Bosnia_and_Herzegovina': 'Europe',
    'Bulgaria': 'Europe',
    'Croatia': 'Europe',
    'Cyprus': 'Europe',
    'Czech_Republic': 'Europe',
    'Denmark': 'Europe',
    'Estonia': 'Europe',
    'Finland': 'Europe',
    'France': 'Europe',
    'Germany': 'Europe',
    'Greece': 'Europe',
    'Hungary': 'Europe',
    'Iceland': 'Europe',
    'Ireland': 'Europe',
    'Italy': 'Europe',
    'Kosovo': 'Europe',
    'Latvia': 'Europe',
    'Liechtenstein': 'Europe',
    'Lithuania': 'Europe',
    'Luxembourg': 'Europe',
    'Macedonia': 'Europe',
    'Malta': 'Europe',
    'Moldova': 'Europe',
    'Monaco': 'Europe',
    'Montenegro': 'Europe',
    'Netherlands': 'Europe',
    'Norway': 'Europe',
    'Poland': 'Europe',
    'Portugal': 'Europe',
    'Romania': 'Europe',
    'Russia': 'Europe',
    'San_Marino': 'Europe',
    'Serbia': 'Europe',
    'Slovakia': 'Europe',
    'Slovenia': 'Europe',
    'Spain': 'Europe',
    'Sweden': 'Europe',
    'Switzerland': 'Europe',
    'Ukraine': 'Europe',
    'United_Kingdom': 'Europe',
    'Vatican_City': 'Europe',
    
    # North America
    'Antigua_and_Barbuda': 'North_America',
    'Bahamas': 'North_America',
    'Barbados': 'North_America',
    'Belize': 'North_America',
    'Canada': 'North_America',
    'Costa_Rica': 'North_America',
    'Cuba': 'North_America',
    'Dominica': 'North_America',
    'Dominican_Republic': 'North_America',
    'El_Salvador': 'North_America',
    'Grenada': 'North_America',
    'Guatemala': 'North_America',
    'Haiti': 'North_America',
    'Honduras': 'North_America',
    'Jamaica': 'North_America',
    'Mexico': 'North_America',
    'Nicaragua': 'North_America',
    'Panama': 'North_America',
    'Saint_Kitts_and_Nevis': 'North_America',
    'Saint_Lucia': 'North_America',
    'Saint_Vincent_and_the_Grenadines': 'North_America',
    'Trinidad_and_Tobago': 'North_America',
    'United_States': 'North_America',
    'USA': 'North_America',
    
    # South America
    'Argentina': 'South_America',
    'Bolivia': 'South_America',
    'Brazil': 'South_America',
    'Chile': 'South_America',
    'Colombia': 'South_America',
    'Ecuador': 'South_America',
    'Guyana': 'South_America',
    'Paraguay': 'South_America',
    'Peru': 'South_America',
    'Suriname': 'South_America',
    'Uruguay': 'South_America',
    'Venezuela': 'South_America',
    
    # Oceania
    'American_Samoa': 'Oceania',
    'Australia': 'Oceania',
    'Fiji': 'Oceania',
    'Kiribati': 'Oceania',
    'Marshall_Islands': 'Oceania',
    'Micronesia': 'Oceania',
    'Nauru': 'Oceania',
    'New_Caledonia': 'Oceania',
    'New_Zealand': 'Oceania',
    'Palau': 'Oceania',
    'Papua_New_Guinea': 'Oceania',
    'Samoa': 'Oceania',
    'Solomon_Islands': 'Oceania',
    'Tonga': 'Oceania',
    'Tuvalu': 'Oceania',
    'Vanuatu': 'Oceania',
    
    # Special territories/regions
    'Hong_Kong': 'Asia',
    'Gibraltar': 'Europe',
    'Viet_Nam': 'Asia',
}

def get_continent(country):
    """Get continent for a country, return Unknown if not found."""
    return COUNTRY_TO_CONTINENT.get(country, 'Unknown')

def replace_country_with_continent(input_fasta, output_fasta):
    """
    Replace country with continent in FASTA headers.
    Header format: >accession/date/country/length/genotype
    """
    
    total_sequences = 0
    converted_sequences = 0
    unknown_countries = set()
    continent_counts = {}
    
    modified_records = []
    
    print("Processing sequences and replacing countries with continents...\n")
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        total_sequences += 1
        
        # Parse header
        header_parts = record.id.split('/')
        
        if len(header_parts) >= 5:
            accession = header_parts[0]
            date = header_parts[1]
            country = header_parts[2]
            length = header_parts[3]
            genotype = header_parts[4]
            
            # Get continent
            continent = get_continent(country)
            
            if continent == 'Unknown':
                unknown_countries.add(country)
            else:
                converted_sequences += 1
            
            # Count continents
            continent_counts[continent] = continent_counts.get(continent, 0) + 1
            
            # Create new header
            new_id = f"{accession}/{date}/{continent}/{length}/{genotype}"
            record.id = new_id
            record.description = ""  # Clear description to avoid duplication
            
            modified_records.append(record)
        else:
            # Keep original if format doesn't match
            modified_records.append(record)
        
        # Progress indicator
        if total_sequences % 1000 == 0:
            print(f"Processed {total_sequences} sequences...")
    
    # Write output
    with open(output_fasta, 'w') as fout:
        SeqIO.write(modified_records, fout, "fasta")
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"CONVERSION SUMMARY")
    print(f"{'='*70}")
    print(f"Total sequences processed: {total_sequences}")
    print(f"Successfully converted: {converted_sequences}")
    print(f"Unknown countries: {len(unknown_countries)}")
    
    print(f"\n{'='*70}")
    print(f"CONTINENT DISTRIBUTION")
    print(f"{'='*70}")
    for continent in sorted(continent_counts.keys()):
        count = continent_counts[continent]
        percentage = count / total_sequences * 100
        print(f"{continent}: {count:6d} sequences ({percentage:5.2f}%)")
    
    if unknown_countries:
        print(f"\n{'='*70}")
        print(f"UNKNOWN COUNTRIES (need to be added to mapping)")
        print(f"{'='*70}")
        for country in sorted(unknown_countries):
            print(f"  {country}")
    
    print(f"\n{'='*70}")
    print(f"Output written to: {output_fasta}")
    print(f"{'='*70}")

def main():
    input_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Aligned_N450_final.fasta"
    output_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Aligned_N450_final_continents.fasta"
    
    replace_country_with_continent(input_fasta, output_fasta)

if __name__ == "__main__":
    main()
